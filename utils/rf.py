from gnomad_hail import *
import pyspark.sql
from pyspark.ml.feature import *
from pyspark.ml.classification import *
from pyspark.ml import *
from pyspark.sql import Row
import hail as hl
import pandas as pd


def run_rf_test(
        mt: hl.MatrixTable,
        output: str = '/tmp'
) -> Tuple[pyspark.ml.PipelineModel, hl.MatrixTable]:
    """
    Runs a dummy test RF on a given MT:
    1. Creates row annotations and labels to run model on
    2. Trains a RF pipeline model (including median imputation of missing values in created annotations)
    3. Saves the RF pipeline model
    4. Applies the model to the MT and prints features importance

    :param MatrixTable mt: Input MT
    :param str output: Output files prefix to save the RF model
    :return: RF model and MatrixTable after applying RF model
    :rtype: (PipelineModel, MatrixTable)
    """

    mt = mt.annotate_rows(
                          feature1 = hl.rand_bool(0.1),
                          feature2 = hl.rand_norm(0.0, 1.0),
                          feature3 = hl.or_missing(hl.rand_bool(0.5), hl.rand_norm(0.0, 1.0)))

    mt = mt.annotate_rows(label=hl.cond(mt['feature1'] & (mt['feature2'] > 0), "TP", "FP"))
    ht = mt.rows()

    def f3stats(ht):
        return ht.aggregate(hl.struct(
            n=hl.agg.count_where(hl.is_defined(ht['feature3'])),
            med=hl.median(hl.agg.collect(ht['feature3']))
        )
        )

    f3_before_imputation = f3stats(ht)
    logger.info('Feature3 defined values before imputation: {}'.format(f3_before_imputation.n))
    logger.info('Feature3 median: {}'.format(f3_before_imputation.med))

    features_to_impute = ['feature3']
    quantiles = get_columns_quantiles(ht, features_to_impute, [0.5])
    quantiles = {k: v[0] for k, v in quantiles.items()}

    logger.info('Features median:\n{}'.format(f'{k}: {v}\n' for k,v in quantiles.items()))
    ht = ht.annotate(
        **{f: hl.or_else(ht[f], quantiles[f]) for f in features_to_impute}
    )
    ht = ht.annotate_globals(medians= quantiles)

    f3_after_imputation = f3stats(ht)
    logger.info('Feature3 defined values after imputation: {}'.format(f3_after_imputation.n))
    logger.info('Feature3 median: {}'.format(f3_after_imputation.med))

    ht = ht.select('label', 'feature1', 'feature2', 'feature3')

    label = 'label'
    features = ['feature1', 'feature2', 'feature3']

    rf_model = train_rf(ht, features, label)
    save_model(rf_model, out_path=output + '/rf.model', overwrite=True)
    rf_model = load_model(output + '/rf.model')

    return rf_model, apply_rf_model(ht, rf_model, features, label)

def get_columns_quantiles(
        ht: hl.Table,
        columns: List[str],
        quantiles: List[float],
        relative_error: int = 0.001
) -> Dict[str, List[float]]:
    """
    Computes approximate quantiles of specified numeric columns from non-missing values.
    Non-numeric fields are ignored.
    This function returns a Dict of column name -> list of quantiles in the same order specified.
    If a column only has NAs, None is returned.

    :param Table ht: input HT
    :param list of str features: list of features to impute. If none given, all numerical features with missing data are imputed
    :param float relative_error: The relative error on the quantile approximation
    :return: Dict of column -> quantiles
    :rtype: dict of str -> list of float
    """

    df = ht.select(*columns).to_spark()

    res = {}
    for f in columns:
        logger.info("Computing median for column: {}".format(f))
        col_no_na = df.select(f).dropna()
        if col_no_na.first() is not None:
            res[f] = col_no_na.approxQuantile(str(f), quantiles, relative_error)
        else:
            res[f] = None

    return res


def ht_to_rf_df(
        ht: hl.Table,
        features: List[str],
        label: str,
        index: str = None
) -> pyspark.sql.DataFrame:
    """

    Creates a Spark dataframe ready for RF from a HT.
    Rows with any missing features are dropped.
    Missing labels are replaced with 'NA'
    Note: Only basic types are supported!

    :param Table ht: Input HT
    :param list of str features: Features that will be used for RF
    :param str label: Label column that will be predicted by RF
    :param str index: Optional index column to keep (E.g. for joining results back at later stage)
    :return: Spark Dataframe
    :rtype: DataFrame
    """

    cols_to_keep = features + [label]
    if index:
        cols_to_keep.append(index)

    df = ht.select(*cols_to_keep).to_spark()
    df = df.dropna(subset=features).fillna('NA', subset=label)

    return df


def get_features_importance(
        rf_pipeline: pyspark.ml.PipelineModel,
        rf_index: int = -2,
        assembler_index: int =-3
) -> Dict[str, float]:
    """
    Extract the features importance from a Pipeline model containing a RandomForestClassifier stage.

    :param PipelineModel rf_pipeline: Input pipeline
    :param int rf_index: index of the RandomForestClassifier stage
    :param int assembler_index: index of the VectorAssembler stage
    :return: feature importance for each feature in the RF model
    :rtype: dict of str: float
    """

    feature_names = [x[:-len("_indexed")] if x.endswith("_indexed") else x for x in
                     rf_pipeline.stages[assembler_index].getInputCols()]
    feature_importance = {fromSSQL(new_name): importance for
                          (new_name, importance) in zip(feature_names, rf_pipeline.stages[rf_index].featureImportances)}
    return feature_importance


def get_labels(
        rf_pipeline: pyspark.ml.PipelineModel
) -> List[str]:
    """
    Returns the labels from the StringIndexer stage at index 0 from an RF pipeline model

    :param PipelineModel rf_pipeline: Input pipeline
    :return: labels
    :rtype: list of str
    """
    return rf_pipeline.stages[0].labels


def test_model(
        ht: hl.Table,
        rf_model: pyspark.ml.PipelineModel,
        features: List[str],
        label: str,
        prediction_col_name: str = 'rf_prediction'
) -> List(hl.tstruct):
    """
    A wrapper to test a model on a set of examples with known labels:
    1) Runs the model on the data
    2) Prints confusion matrix and accuracy
    3) Returns confusion matrix as a list of struct

    :param Table ht: Input table
    :param PipelineModel rf_model: RF Model
    :param list of str features: Columns containing features that were used in the model
    :param str label: Column containing label to be predicted
    :param str prediction_col_name: Where to store the prediction
    :return: A list containing structs with {label, prediction, n}
    :rtype: list of Struct
    """

    ht = apply_rf_model(ht.filter(hl.is_defined(ht[label])),
                        rf_model,
                        features,
                        label,
                        prediction_col_name=prediction_col_name)

    test_results = ht.group_by(ht[prediction_col_name], ht[label]).aggregate(n=hl.agg.count()).collect()

    #Print results
    df = pd.DataFrame(test_results)
    df = df.pivot(index=label, columns=prediction_col_name, values='n')
    logger.info("Testing results:\n{}".format(pformat(df)))
    logger.info("Accuracy: {}".format(
        sum([x.n for x in test_results if x[label] == x['rf_prediction']]) /
        sum([x.n for x in test_results])
    ))

    return test_results



def apply_rf_model(
        ht: hl.Table,
        rf_model: pyspark.ml.PipelineModel,
        features: List[str],
        label: str,
        probability_col_name: str = 'rf_probability',
        prediction_col_name: str = 'rf_prediction'
) -> hl.Table:
    """
    Applies a Random Forest (RF) pipeline model to a Table and annotate the RF probabilities and predictions.

    :param MatrixTable ht: Input HT
    :param PipelineModel rf_model: Random Forest pipeline model
    :param list of str features: List of feature columns in the pipeline. !Should match the model list of features!
    :param str label: Column containing the labels. !Should match the model labels!
    :param str probability_col_name: Name of the column that will store the RF probabilities
    :param str prediction_col_name: Name of the column that will store the RF predictions
    :return: Table with RF columns
    :rtype: Table
    """

    logger.info("Applying RF model.")

    index_name = 'rf_idx'
    while index_name in ht.row:
        index_name += '_tmp'
    ht = ht.add_index(name=index_name)

    ht_keys = ht.key
    ht = ht.key_by(index_name)

    df = ht_to_rf_df(ht, features, label, index_name)

    rf_df = rf_model.transform(df)

    rf_ht = hl.Table.from_spark(
        rf_df.rdd.map(
            lambda row:
            Row(idx=row[index_name],
                probability=row["probability"].toArray().tolist(),
                predictedLabel=row["predictedLabel"])
        ).toDF()
    ).persist()

    rf_ht = rf_ht.key_by('idx')

    ht = ht.annotate(
        **{
            probability_col_name: {label: rf_ht[ht[index_name]]["probability"][i] for i, label in enumerate(get_labels(rf_model))},
            prediction_col_name: rf_ht[ht[index_name]]["predictedLabel"]
        }
    )

    ht = ht.key_by(*ht_keys)
    ht = ht.drop(index_name)

    return ht


def save_model(
        rf_pipeline: pyspark.ml.PipelineModel,
        out_path: str,
        overwrite: bool = False
) -> None:
    """
    Saves a Random Forest pipeline model.

    :param PipelineModel rf_pipeline: Pipeline to save
    :param str out_path: Output path
    :param bool overwrite: If set, will overwrite existing file(s) at output location
    :return: Nothing
    :rtype: NoneType
    """
    logger.info("Saving model to %s" % out_path)
    if overwrite:
        rf_pipeline.write().overwrite().save(out_path)
    else:
        rf_pipeline.save(out_path)


def load_model(
        input_path: str
) -> pyspark.ml.PipelineModel:
    """
    Loads a Random Forest pipeline model.

    :param str input_path: Location of model to load
    :return: Random Forest pipeline model
    :rtype: PipelineModel
    """
    logger.info("Loading model from {}".format(input_path))
    return pyspark.ml.PipelineModel.load(input_path)


def train_rf(
        ht: hl.Table,
        features: List[str],
        label: str,
        num_trees: int = 500,
        max_depth: int = 5
) -> pyspark.ml.PipelineModel:
    """
    Trains a Random Forest (RF) pipeline model.

    :param Table ht: Input HT
    :param list of str features: List of columns to be used as features
    :param str label: Column containing the label to predict
    :param int num_trees: Number of trees to use
    :param int max_depth: Maximum tree depth
    :return: Random Forest pipeline model
    :rtype: PipelineModel
    """
    label_name = list(ht.select(label).row)[0]
    feature_names = ht.select(*features).columns

    logger.info("Training RF model using:\n"
                "features: {}\n"
                "labels: {}\n"
                "num_trees: {}\n"
                "max_depth: {}".format(",".join(feature_names),
                                   label_name, num_trees, max_depth))

    df = ht_to_rf_df(ht, features, label)

    label_indexer = StringIndexer(inputCol=label_name, outputCol=label_name + "_indexed").fit(df)
    labels = label_indexer.labels
    logger.info("Found labels: {}".format(labels))

    string_features = [x[0] for x in df.dtypes if x[0] != label_name and x[1] == 'string']
    if string_features:
        logger.info("Indexing string features: {}".format(",".join(string_features)))
    string_features_indexers = [StringIndexer(inputCol=x, outputCol=x + "_indexed").fit(df)
                                for x in string_features]

    assembler = VectorAssembler(inputCols=[x[0] + "_indexed" if x[1] == 'string' else x[0]
                                           for x in df.dtypes if x[0] != label_name ],
                                outputCol="features")

    rf = RandomForestClassifier(labelCol=label_name + "_indexed", featuresCol="features",
                                maxDepth=max_depth, numTrees=num_trees)

    label_converter = IndexToString(inputCol='prediction', outputCol='predictedLabel', labels=labels)

    pipeline = Pipeline(stages=[label_indexer] + string_features_indexers +
                               [assembler, rf, label_converter])

    #Train model
    logger.info("Training RF model")
    rf_model = pipeline.fit(df)

    feature_importance = get_features_importance(rf_model)

    logger.info(
        "RF features importance:\n{}".format(
            "\n".join(["{}: {}".format(f, i) for f, i in feature_importance.items()])
        )
    )

    return rf_model
