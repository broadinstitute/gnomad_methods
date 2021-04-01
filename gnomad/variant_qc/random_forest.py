import json
import logging
import pprint
from pprint import pformat
from typing import Dict, List, Optional, Tuple

import hail as hl
import pandas as pd
import pyspark.sql
from pyspark.ml import Pipeline
from pyspark.ml.classification import RandomForestClassifier
from pyspark.ml.feature import IndexToString, StringIndexer, VectorAssembler
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf  # pylint: disable=no-name-in-module
from pyspark.sql.types import ArrayType, DoubleType

from gnomad.utils.file_utils import file_exists

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def run_rf_test(
    mt: hl.MatrixTable, output: str = "/tmp"
) -> Tuple[pyspark.ml.PipelineModel, hl.Table]:
    """
    Runs a dummy test RF on a given MT.

    1. Creates row annotations and labels to run model on
    2. Trains a RF pipeline model (including median imputation of missing values in created annotations)
    3. Saves the RF pipeline model
    4. Applies the model to the MT and prints features importance

    :param mt: Input MT
    :param output: Output files prefix to save the RF model
    :return: RF model and MatrixTable after applying RF model
    """

    mt = mt.annotate_rows(
        feature1=hl.rand_bool(0.1),
        feature2=hl.rand_norm(0.0, 1.0),
        feature3=hl.or_missing(hl.rand_bool(0.5), hl.rand_norm(0.0, 1.0)),
    )

    mt = mt.annotate_rows(
        label=hl.cond(mt["feature1"] & (mt["feature2"] > 0), "TP", "FP")
    )
    ht = mt.rows()

    def f3stats(ht):
        return ht.aggregate(
            hl.struct(
                n=hl.agg.count_where(hl.is_defined(ht["feature3"])),
                med=hl.median(hl.agg.collect(ht["feature3"])),
            )
        )

    f3_before_imputation = f3stats(ht)
    logger.info(
        "Feature3 defined values before imputation: {}".format(f3_before_imputation.n)
    )
    logger.info("Feature3 median: {}".format(f3_before_imputation.med))

    features_to_impute = ["feature3"]
    quantiles = get_columns_quantiles(ht, features_to_impute, [0.5])
    quantiles = {k: v[0] for k, v in quantiles.items()}

    logger.info(
        "Features median:\n{}".format(f"{k}: {v}\n" for k, v in quantiles.items())
    )
    ht = ht.annotate(**{f: hl.or_else(ht[f], quantiles[f]) for f in features_to_impute})
    ht = ht.annotate_globals(medians=quantiles)

    f3_after_imputation = f3stats(ht)
    logger.info(
        "Feature3 defined values after imputation: {}".format(f3_after_imputation.n)
    )
    logger.info("Feature3 median: {}".format(f3_after_imputation.med))

    ht = ht.select("label", "feature1", "feature2", "feature3")

    label = "label"
    features = ["feature1", "feature2", "feature3"]

    rf_model = train_rf(ht, features, label)
    save_model(rf_model, out_path=output + "/rf.model", overwrite=True)
    rf_model = load_model(output + "/rf.model")

    return rf_model, apply_rf_model(ht, rf_model, features, label)


def check_ht_fields_for_spark(ht: hl.Table, fields: List[str]) -> None:
    """
    Checks specified fields of a hail table for Spark DataFrame conversion (type and name)

    :param ht: input Table
    :param fields: Fields to test
    :return: None
    """

    allowed_types = [
        hl.tfloat,
        hl.tfloat32,
        hl.tfloat64,
        hl.tint,
        hl.tint32,
        hl.tint64,
        hl.tstr,
        hl.tbool,
    ]

    bad_field_names = [c for c in fields if "." in c]

    bad_types = [
        c[0]
        for c in ht.key_by().select(*fields).row.items()
        if c[1].dtype not in allowed_types
    ]

    if bad_field_names or bad_types:
        raise ValueError(
            "Only basic type fields can be converted from Hail to Spark. In addition, `.` are not allowed in field names in Spark.\n"
            + "Offending fields (non basic type): {}".format(bad_types)
            + "Offending fields (bad field name): {}\n".format(
                ",".join(bad_field_names)
            )
        )

    return


def get_columns_quantiles(
    ht: hl.Table, fields: List[str], quantiles: List[float], relative_error: int = 0.001
) -> Dict[str, List[float]]:
    """
    Computes approximate quantiles of specified numeric fields from non-missing values. Non-numeric fields are ignored.

    This function returns a Dict of column name -> list of quantiles in the same order specified.
    If a column only has NAs, None is returned.

    :param ht: input HT
    :param fields: list of features to impute. If none given, all numerical features with missing data are imputed
    :param quantiles: list of quantiles to return (e.g. [0.5] would return the median)
    :param relative_error: The relative error on the quantile approximation
    :return: Dict of column -> quantiles
    """

    check_ht_fields_for_spark(ht, fields)

    df = ht.key_by().select(*fields).to_spark()

    res = {}
    for f in fields:
        logger.info("Computing median for column: {}".format(f))
        col_no_na = df.select(f).dropna()
        if col_no_na.first() is not None:
            res[f] = col_no_na.approxQuantile(str(f), quantiles, relative_error)
        else:
            res[f] = None

    return res


def median_impute_features(
    ht: hl.Table, strata: Optional[Dict[str, hl.expr.Expression]] = None
) -> hl.Table:
    """
    Numerical features in the Table are median-imputed by Hail's `approx_median`.

    If a `strata` dict is given, imputation is done based on the median of of each stratification.

    The annotations that are added to the Table are
        - feature_imputed - A row annotation indicating if each numerical feature was imputed or not.
        - features_median - A global annotation containing the median of the numerical features. If `strata` is given,
          this struct will also be broken down by the given strata.
        - variants_by_strata - An additional global annotation with the variant counts by strata that will only be
          added if imputing by a given `strata`.

    :param ht: Table containing all samples and features for median imputation.
    :param strata: Whether to impute features median by specific strata (default False).
    :return: Feature Table imputed using approximate median values.
    """

    logger.info("Computing feature medians for imputation of missing numeric values")
    numerical_features = [
        k for k, v in ht.row.dtype.items() if v == hl.tint or v == hl.tfloat
    ]

    median_agg_expr = hl.struct(
        **{feature: hl.agg.approx_median(ht[feature]) for feature in numerical_features}
    )

    if strata:
        ht = ht.annotate_globals(
            feature_medians=ht.aggregate(
                hl.agg.group_by(hl.tuple([ht[x] for x in strata]), median_agg_expr),
                _localize=False,
            ),
            variants_by_strata=ht.aggregate(
                hl.agg.counter(hl.tuple([ht[x] for x in strata])), _localize=False
            ),
        )
        feature_median_expr = ht.feature_medians[hl.tuple([ht[x] for x in strata])]
        logger.info(
            "Variant count by strata:\n{}".format(
                "\n".join(
                    [
                        "{}: {}".format(k, v)
                        for k, v in hl.eval(ht.variants_by_strata).items()
                    ]
                )
            )
        )

    else:
        ht = ht.annotate_globals(
            feature_medians=ht.aggregate(median_agg_expr, _localize=False)
        )
        feature_median_expr = ht.feature_medians

    ht = ht.annotate(
        **{f: hl.or_else(ht[f], feature_median_expr[f]) for f in numerical_features},
        feature_imputed=hl.struct(
            **{f: hl.is_missing(ht[f]) for f in numerical_features}
        ),
    )

    return ht


def ht_to_rf_df(
    ht: hl.Table, features: List[str], label: str, index: str = None
) -> pyspark.sql.DataFrame:
    """
    Creates a Spark dataframe ready for RF from a HT.
    Rows with any missing features are dropped.
    Missing labels are replaced with 'NA'

    .. note::

        Only basic types are supported!

    :param ht: Input HT
    :param features: Features that will be used for RF
    :param label: Label column that will be predicted by RF
    :param index: Optional index column to keep (E.g. for joining results back at later stage)
    :return: Spark Dataframe
    """

    cols_to_keep = features + [label]
    if index:
        cols_to_keep.append(index)

    df = ht.key_by().select(*cols_to_keep).to_spark()
    df = df.dropna(subset=features).fillna("NA", subset=label)

    return df


def get_features_importance(
    rf_pipeline: pyspark.ml.PipelineModel, rf_index: int = -2, assembler_index: int = -3
) -> Dict[str, float]:
    """
    Extract the features importance from a Pipeline model containing a RandomForestClassifier stage.

    :param rf_pipeline: Input pipeline
    :param rf_index: index of the RandomForestClassifier stage
    :param assembler_index: index of the VectorAssembler stage
    :return: feature importance for each feature in the RF model
    """

    feature_names = [
        x[: -len("_indexed")] if x.endswith("_indexed") else x
        for x in rf_pipeline.stages[assembler_index].getInputCols()
    ]

    return dict(zip(feature_names, rf_pipeline.stages[rf_index].featureImportances))


def get_labels(rf_pipeline: pyspark.ml.PipelineModel) -> List[str]:
    """
    Returns the labels from the StringIndexer stage at index 0 from an RF pipeline model

    :param rf_pipeline: Input pipeline
    :return: labels
    """
    return rf_pipeline.stages[0].labels


def test_model(
    ht: hl.Table,
    rf_model: pyspark.ml.PipelineModel,
    features: List[str],
    label: str,
    prediction_col_name: str = "rf_prediction",
) -> List[hl.tstruct]:
    """
    A wrapper to test a model on a set of examples with known labels.

    1) Runs the model on the data
    2) Prints confusion matrix and accuracy
    3) Returns confusion matrix as a list of struct

    :param ht: Input table
    :param rf_model: RF Model
    :param features: Columns containing features that were used in the model
    :param label: Column containing label to be predicted
    :param prediction_col_name: Where to store the prediction
    :return: A list containing structs with {label, prediction, n}
    """

    ht = apply_rf_model(
        ht.filter(hl.is_defined(ht[label])),
        rf_model,
        features,
        label,
        prediction_col_name=prediction_col_name,
    )

    test_results = (
        ht.group_by(ht[prediction_col_name], ht[label])
        .aggregate(n=hl.agg.count())
        .collect()
    )

    # Print results
    df = pd.DataFrame(test_results)
    df = df.pivot(index=label, columns=prediction_col_name, values="n")
    logger.info("Testing results:\n{}".format(pprint.pformat(df)))
    logger.info(
        "Accuracy: {}".format(
            sum([x.n for x in test_results if x[label] == x[prediction_col_name]])
            / sum([x.n for x in test_results])
        )
    )

    return test_results


def apply_rf_model(
    ht: hl.Table,
    rf_model: pyspark.ml.PipelineModel,
    features: List[str],
    label: str,
    probability_col_name: str = "rf_probability",
    prediction_col_name: str = "rf_prediction",
) -> hl.Table:
    """
    Applies a Random Forest (RF) pipeline model to a Table and annotate the RF probabilities and predictions.

    :param ht: Input HT
    :param rf_model: Random Forest pipeline model
    :param features: List of feature columns in the pipeline. !Should match the model list of features!
    :param label: Column containing the labels. !Should match the model labels!
    :param probability_col_name: Name of the column that will store the RF probabilities
    :param prediction_col_name: Name of the column that will store the RF predictions
    :return: Table with RF columns
    """

    logger.info("Applying RF model.")

    check_ht_fields_for_spark(ht, features + [label])

    index_name = "rf_idx"
    while index_name in ht.row:
        index_name += "_tmp"
    ht = ht.add_index(name=index_name)

    ht_keys = ht.key
    ht = ht.key_by(index_name)

    df = ht_to_rf_df(ht, features, label, index_name)

    rf_df = rf_model.transform(df)

    def to_array(col):
        def to_array_(v):
            return v.toArray().tolist()

        return udf(to_array_, ArrayType(DoubleType()))(col)

    # Note: SparkSession is needed to write DF to disk before converting to HT;
    # hail currently fails without intermediate write
    spark = SparkSession.builder.getOrCreate()
    rf_df.withColumn("probability", to_array(col("probability"))).select(
        [index_name, "probability", "predictedLabel"]
    ).write.mode("overwrite").save("rf_probs.parquet")
    rf_df = spark.read.format("parquet").load("rf_probs.parquet")
    rf_ht = hl.Table.from_spark(rf_df)
    rf_ht = rf_ht.checkpoint("/tmp/rf_raw_pred.ht", overwrite=True)
    rf_ht = rf_ht.key_by(index_name)

    ht = ht.annotate(
        **{
            probability_col_name: {
                label: rf_ht[ht[index_name]]["probability"][i]
                for i, label in enumerate(get_labels(rf_model))
            },
            prediction_col_name: rf_ht[ht[index_name]]["predictedLabel"],
        }
    )

    ht = ht.key_by(*ht_keys)
    ht = ht.drop(index_name)

    return ht


def save_model(
    rf_pipeline: pyspark.ml.PipelineModel, out_path: str, overwrite: bool = False
) -> None:
    """
    Saves a Random Forest pipeline model.

    :param rf_pipeline: Pipeline to save
    :param out_path: Output path
    :param overwrite: If set, will overwrite existing file(s) at output location
    :return: Nothing
    """
    logger.info("Saving model to %s" % out_path)
    if overwrite:
        rf_pipeline.write().overwrite().save(out_path)
    else:
        rf_pipeline.save(out_path)


def load_model(input_path: str) -> pyspark.ml.PipelineModel:
    """
    Loads a Random Forest pipeline model.

    :param input_path: Location of model to load
    :return: Random Forest pipeline model
    """
    logger.info("Loading model from {}".format(input_path))
    return pyspark.ml.PipelineModel.load(input_path)


def train_rf(
    ht: hl.Table,
    features: List[str],
    label: str,
    num_trees: int = 500,
    max_depth: int = 5,
) -> pyspark.ml.PipelineModel:
    """
    Trains a Random Forest (RF) pipeline model.

    :param ht: Input HT
    :param features: List of columns to be used as features
    :param label: Column containing the label to predict
    :param num_trees: Number of trees to use
    :param max_depth: Maximum tree depth
    :return: Random Forest pipeline model
    """

    logger.info(
        "Training RF model using:\n"
        "features: {}\n"
        "labels: {}\n"
        "num_trees: {}\n"
        "max_depth: {}".format(",".join(features), label, num_trees, max_depth)
    )

    check_ht_fields_for_spark(ht, features + [label])

    df = ht_to_rf_df(ht, features, label)

    label_indexer = (
        StringIndexer(inputCol=label, outputCol=label + "_indexed")
        .setHandleInvalid("keep")
        .fit(df)
    )
    labels = label_indexer.labels
    logger.info("Found labels: {}".format(labels))

    string_features = [x[0] for x in df.dtypes if x[0] != label and x[1] == "string"]
    if string_features:
        logger.info("Indexing string features: {}".format(",".join(string_features)))
    string_features_indexers = [
        StringIndexer(inputCol=x, outputCol=x + "_indexed")
        .setHandleInvalid("keep")
        .fit(df)
        for x in string_features
    ]

    assembler = VectorAssembler(
        inputCols=[
            x[0] + "_indexed" if x[1] == "string" else x[0]
            for x in df.dtypes
            if x[0] != label
        ],
        outputCol="features",
    )

    rf = RandomForestClassifier(
        labelCol=label + "_indexed",
        featuresCol="features",
        maxDepth=max_depth,
        numTrees=num_trees,
    )

    label_converter = IndexToString(
        inputCol="prediction", outputCol="predictedLabel", labels=labels
    )

    pipeline = Pipeline(
        stages=[label_indexer]
        + string_features_indexers
        + [assembler, rf, label_converter]
    )

    # Train model
    logger.info("Training RF model")
    rf_model = pipeline.fit(df)

    feature_importance = get_features_importance(rf_model)

    logger.info(
        "RF features importance:\n{}".format(
            "\n".join(["{}: {}".format(f, i) for f, i in feature_importance.items()])
        )
    )

    return rf_model


def get_rf_runs(rf_json_fp: str) -> Dict:
    """
    Loads RF run data from JSON file.

    :param rf_json_fp: File path to rf json file.
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if file_exists(rf_json_fp):
        with hl.hadoop_open(rf_json_fp) as f:
            return json.load(f)
    else:
        logger.warning(
            f"File {rf_json_fp} could not be found. Returning empty RF run hash dict."
        )
        return {}


def get_run_data(
    input_args: Dict[str, bool],
    test_intervals: List[str],
    features_importance: Dict[str, float],
    test_results: List[hl.tstruct],
) -> Dict:
    """
    Creates a Dict containing information about the RF input arguments and feature importance

    :param Dict of bool keyed by str input_args: Dictionary of model input arguments
    :param List of str test_intervals: Intervals withheld from training to be used in testing
    :param Dict of float keyed by str features_importance: Feature importance returned by the RF
    :param List of struct test_results: Accuracy results from applying RF model to the test intervals
    :return: Dict of RF information
    """
    run_data = {
        "input_args": input_args,
        "features_importance": features_importance,
        "test_intervals": test_intervals,
    }

    if test_results is not None:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            # Note: values[0] is the TP/FP label and values[1] is the prediction
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data["test_results"] = [dict(x) for x in test_results]
        run_data["test_accuracy"] = tps / total

    return run_data


def pretty_print_runs(
    runs: Dict, label_col: str = "rf_label", prediction_col_name: str = "rf_prediction"
) -> None:
    """
    Prints the information for the RF runs loaded from the json file storing the RF run hashes -> info

    :param runs: Dictionary containing JSON input loaded from RF run file
    :param label_col: Name of the RF label column
    :param prediction_col_name: Name of the RF prediction column
    :return: Nothing -- only prints information
    """

    for run_hash, run_data in runs.items():
        print(f"\n=== {run_hash} ===")
        testing_results = (
            run_data.pop("test_results") if "test_results" in run_data else None
        )
        # print(testing_results)
        print(json.dumps(run_data, sort_keys=True, indent=4, separators=(",", ": ")))
        if testing_results is not None:
            # Print results
            res_pd = pd.DataFrame(testing_results)
            res_pd = res_pd.pivot(
                index=label_col, columns=prediction_col_name, values="n"
            )
            logger.info("Testing results:\n{}".format(pformat(res_pd)))
