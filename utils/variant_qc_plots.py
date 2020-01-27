from bokeh.palettes import d3  # pylint: disable=no-name-in-module

from gnomad_hail import *
from gnomad_qc.resources.variant_qc import score_ranking_path, rf_path, binned_concordance_path  # FIXME: Need to update so there is no reference to gnomad_qc
from gnomad_hail.utils.plotting import *

qc_plots_settings = {
    'mean_point_size': 4.0,
    'min_point_size': 1.0,
    'max_point_size': 16.0,
    'label_text_font_size': "14pt",
    'title.text_font_size': "16pt",
    'subtitle.text_font_size': "14pt",
    'axis.axis_label_text_font_size': "16pt",
    'axis.axis_label_text_font_style': "normal",
    'axis.major_label_text_font_size': "14pt"
}


def get_point_size_col(data: pd.Series, size_prop: str) -> pd.Series:
    """
    Given a data Series, returns the corresponding point size either:
    - Constant to qc_plots_settings['mean_point_size'] if `size_prop` is None
    - Radius proportional to data, if `size_prop` is 'radius'
    - Area proportional to data, if `size_prop` is 'area'

    Mean, min and max point  sizes are extracted from qc_plots_settings

    :param Series data: Input data series
    :param str size_prop: One of None, 'radius' or 'area'
    :return: Series with corresponding point size for each data point
    :rtype: Series
    """
    if size_prop is None:
        return pd.Series(len(data) * [qc_plots_settings['mean_point_size']])
    else:
        mean_data = np.mean(data)
        if size_prop == 'radius':
            return data.apply(lambda x: max(qc_plots_settings['min_point_size'], min(qc_plots_settings['max_point_size'], qc_plots_settings['mean_point_size'] * (x / mean_data))))
        elif size_prop == 'area':
            return data.apply(lambda x: max(qc_plots_settings['min_point_size'], min(qc_plots_settings['max_point_size'], qc_plots_settings['mean_point_size'] * np.pi * (np.sqrt(x / mean_data) / np.pi))))
        else:
            raise ValueError(f"{size_prop} is not a supported value for argument `size_prop`")


def set_plots_defaults(p: Plot) -> None:
    p.legend.label_text_font_size = qc_plots_settings['label_text_font_size']
    p.title.text_font_size = qc_plots_settings['title.text_font_size']
    p.axis.axis_label_text_font_size = qc_plots_settings['axis.axis_label_text_font_size']
    p.axis.axis_label_text_font_style = qc_plots_settings['axis.axis_label_text_font_style']
    p.axis.major_label_text_font_size = qc_plots_settings['axis.major_label_text_font_size']


def get_binned_models_pd(data_type: str, models: Union[Dict[str, str], List[str]], contigs: Set[str] = None) -> pd.DataFrame:
    """
    Creates a single DataFrame with all desired models binned and ready for plotting.

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str models: Models to load. Either a list of the model ids, or a dict with model id -> model name for display
    :param list of str contigs: Contigs to load
    :return: Plot-ready DataFrame
    :rtype: DataFrame
    """

    def aggregate_contig(ht: hl.Table, contigs: Set[str] = None):
        """
        Aggregates all contigs together and computes number for bins accross the contigs.
        """
        if contigs:
            ht = ht.filter(hl.literal(contigs).contains(ht.contig))

        return ht.group_by(
            *[k for k in ht.key if k != 'contig']
        ).aggregate(
            min_score=hl.agg.min(ht.min_score),
            max_score=hl.agg.max(ht.max_score),
            **{x: hl.agg.sum(ht[x]) for x in ht.row_value if x not in ['min_score', 'max_score']}
        )

    if not isinstance(models, dict):
        models = {m: m for m in models}

    hts = [
        aggregate_contig(hl.read_table(score_ranking_path(data_type, model_id, binned=True)), contigs)
            .annotate(model=model_name) for model_id, model_name in models.items()
    ]

    ht = hts[0].union(*hts[1:])
    ht = ht.annotate(n_biallelic=hl.cond(ht.bi_allelic, ht.n, 0))
    return ht.to_pandas()


def plot_metric(df: pd.DataFrame,
                y_name: str,
                cols: List[str],
                y_fun: Callable[[pd.Series], Union[float, int]] = lambda x: x,
                cut: int = None,
                plot_all: bool = True,
                plot_bi_allelics: bool = True,
                plot_singletons: bool = True,
                plot_bi_allelic_singletons: bool = True,
                plot_adj: bool = False,
                colors: Dict[str, str] = None,
                link_cumul_y: bool = True,
                size_prop: str = 'area'
                ) -> Tabs:
    """
    Generic function for generating QC metric plots using a plotting-ready DataFrame (obtained from `get_binned_models_pd`)
    DataFrame needs to have a `rank_id` column, a `bin` column and a `model` column (contains the model name and needs to be added to binned table(s))

    This function generates scatter plots with the metric bin on x-axis and a user-defined function on the y-axis.
    The data for the y-axis function needs to from the columns specified in `cols`. The function is specified with the `y_fun` argument and data columns are access as a list.
    As an example, plotting Transition to transversion ratio is done as follows::

        plot_metric(snvs, 'Ti/Tv', ['n_ti', 'n_tv'], y_fun=lambda x: x[0]/x[1], colors=colors)

    In this command, `x[0]` correspond to the  first column selected (`'n_ti'`)  and `x[1]` to the second (`'n_tv'`).


    This function plots a tab for each of the plot condition(s) selected: all, bi-allelics, bi-allelic singletons.
    Within each tab, each row contains a non-cumulative and a cumulative plot of the bins / values.
    If `plot_adj` is set, then an extra row is added plotting only variants in release samples where AC_ADJ>0. The bin for these sites is computed based on those variants only.

    :param pd.DataFrame df: Input data
    :param str y_name: Name of the metric plotted on the y-axis
    :param list of str cols: Columns used to compute the metric plotted
    :param callable y_fun: Function to apply to the columns to generate the metric
    :param int cut: Where to draw the bin cut
    :param bool plot_all: Whether to plot a tab with all variants
    :param bool plot_bi_allelics: Whether to plot a tab with bi-allelic variants only
    :param bool plot_singletons: Whether to plot a tab with singleton variants only
    :param bool plot_bi_allelic_singletons:  Whether to plot a tab with bi-allelic singleton variants only
    :param bool plot_adj: Whether to plot additional rows with adj variants in release samples only
    :param dict of str -> str colors: Mapping of model name -> color
    :param bool link_cumul_y: If set, y-axes of cumulative and non-cumulative plots are linked
    :param str size_prop: Either 'size' or 'area' can be specified. If either is specified, the points will be sized proportionally to the amount of data in that point.
    :return: Plot
    :rtype: Tabs
    """

    def get_row(df: pd.DataFrame, y_name: str, cols: List[str], y_fun: Callable[[pd.Series], Union[float, int]], titles: List[str], link_cumul_y: bool, cut: int = None) -> Row:
        """
        Generates a single row with two plots: a regular scatter plot and a cumulative one.
        Both plots have bins on the x-axis. The y-axis is computed by applying the function `y_fun` on the columns `cols`.

        Data source is shared between the two plots so that highlighting / selection is linked.
        X-axis is shared between the two plots.
        Y-axus is shared if `link_cumul_y` is `True`

        """

        def get_plot(data_source: ColumnDataSource, y_name: str, y_col_name: str, titles: List[str], data_ranges: Tuple[DataRange1d, DataRange1d], cut: int = None) -> Plot:
            """
            Generates a single scatter plot panel
            """

            p = figure(
                title=titles[0],
                x_axis_label='bin',
                y_axis_label=y_name,
                tools="save,pan,box_zoom,reset,wheel_zoom,box_select,lasso_select,help,hover")
            p.x_range = data_ranges[0]
            p.y_range = data_ranges[1]

            if cut:
                p.add_layout(Span(location=cut, dimension='height', line_color='red', line_dash='dashed'))

            # Add circles layouts one model at a time, so that no default legend is generated.
            # Because data is in the same ColumnDataSource, use a BooleanFilter to plot each model separately
            circles = []
            for model in set(data_source.data['model']):
                view = CDSView(source=data_source, filters=[BooleanFilter([x == model for x in data_source.data['model']])])
                circles.append((model, [p.circle('bin', y_col_name, color='_color', size='_size', source=data_source, view=view)]))

            p.select_one(HoverTool).tooltips = [('model', '@model'),
                                                ('bin', '@bin'),
                                                (y_name, f'@{y_col_name}'),
                                                ('min_score', '@min_score'),
                                                ('max_score', '@max_score'),
                                                ('n_data_points', '@_n')
                                                ] + [(col, f'@{col}') for col in cols]
            set_plots_defaults(p)

            # Add legend above the plot area
            legend = Legend(items=circles, orientation='horizontal', location=(0, 0), click_policy="hide")
            p.add_layout(legend, 'above')

            # Add subtitles if any
            for title in titles[1:]:
                p.add_layout(Title(text=title, text_font_size=qc_plots_settings['subtitle.text_font_size']), 'above')

            return p

        # Compute non-cumulative values by applying `y_fun`
        df['non_cumul'] = df[cols].apply(y_fun, axis=1)

        # Compute cumulative values for each of the data columns
        for col in cols:
            df[f'{col}_cumul'] = df.groupby('model').aggregate(np.cumsum)[col]
        df['cumul'] = df[[f'{col}_cumul' for col in cols]].apply(y_fun, axis=1)

        # Create data ranges that are either shared or distinct depending on the y_cumul parameter
        non_cumul_data_ranges = (DataRange1d(), DataRange1d())
        cumul_data_ranges = non_cumul_data_ranges if link_cumul_y else (non_cumul_data_ranges[0], DataRange1d())
        data_source = ColumnDataSource(df)

        return Row(get_plot(data_source, y_name, 'non_cumul', titles, non_cumul_data_ranges, cut),
                   get_plot(data_source, y_name, 'cumul', [titles[0] + ', cumulative'] + titles[1:], cumul_data_ranges, cut))

    def prepare_pd(df: pd.DataFrame, cols: List[str], colors: Dict[str, str] = {}, size_prop: str = None):
        """
        Groups a pandas DataFrame by model and bin while keeping relevant columns only.
        Adds 3 columns used for plotting:
        1. A _color column column
        2. A _n column containing the number of data points
        3. A _size column containing the size of data points based on the `size_prop` and `qc_plot_settings` parameters
        """
        df = df.groupby(['model', 'bin']).agg({**{col: np.sum for col in cols},
                                               'min_score': np.min, 'max_score': np.max})
        df = df.reset_index()
        df['_color'] = [colors.get(x, 'gray') for x in df['model']]
        df['_n'] = np.sum(df[cols], axis=1)
        df['_size'] = get_point_size_col(df['_n'], size_prop)
        return df

    colors = colors if colors is not None else {}
    tabs = []
    adj_strats = ['', 'adj_'] if plot_adj else ['']

    if plot_all:
        children = []
        for adj in adj_strats:
            titles = [y_name, 'Adj variants (adj rank)' if adj else 'All variants']
            plot_df = prepare_pd(df.loc[df.rank_id == f'{adj}rank'], cols, colors, size_prop)
            if len(plot_df) > 0:
                children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
            else:
                logger.warn('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='All'))

    if plot_bi_allelics:
        children = []
        for adj in adj_strats:
            for biallelic_rank in ['', 'biallelic_']:
                titles = [y_name,'Bi-allelic variants ({} rank)'.format('overall' if not adj and not biallelic_rank else ' '.join([adj[:-1], biallelic_rank[:-1]]).lstrip())]
                plot_df = prepare_pd(df.loc[df.bi_allelic & (df.rank_id == f'{adj}{biallelic_rank}rank')], cols, colors, size_prop)
                if len(plot_df) > 0:
                    children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
                else:
                    logger.warn('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='Bi-allelic'))

    if plot_singletons:
        children = []
        for adj in adj_strats:
            for singleton_rank in ['', 'singleton_']:
                titles = [y_name, 'Singletons ({} rank)'.format('overall' if not adj and not singleton_rank else " ".join([adj[:-1], singleton_rank[:-1]]).lstrip())]
                plot_df = prepare_pd(df.loc[df.singleton & (df.rank_id == f'{adj}{singleton_rank}rank')], cols, colors, size_prop)
                if len(plot_df) > 0:
                    children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
                else:
                    logger.warn('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='Singletons'))

    if plot_bi_allelic_singletons:
        children = []
        for adj in adj_strats:
            for bisingleton_rank in ['', 'biallelic_singleton_']:
                titles = [y_name, 'Bi-allelic singletons ({} rank)'.format('overall' if not adj and not bisingleton_rank else " ".join([adj[:-1], bisingleton_rank[:-1].replace("_", " ")]).lstrip())]
                plot_df = prepare_pd(df.loc[df.bi_allelic & df.singleton & (df.rank_id == f'{adj}{bisingleton_rank}rank')], cols, colors, size_prop)
                if len(plot_df) > 0:
                    children.append(get_row(plot_df, y_name, cols, y_fun, titles, link_cumul_y, cut))
                else:
                    logger.warn('No data found for plot: {}'.format('\t'.join(titles)))

        if children:
            tabs.append(Panel(child=Column(children=children), title='Bi-allelic singletons'))

    return Tabs(tabs=tabs)


def plot_score_distributions(data_type, models: Union[Dict[str, str], List[str]], snv: bool, cut: int, colors: Dict[str, str] = None) -> Tabs:
    """
    Generates plots of model scores distributions:
    One tab per model.
    Within each tab, there is 2x2 grid of plots:
    - One row showing the score distribution across the entire data
    - One row showing the score distribution across the release-samples, adj data only (release_sample_AC_ADJ > 0)
    - One column showing the histogram of the score
    - One column showing the normalized cumulative histogram of the score

    Cutoff is highlighted by a dashed red line

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str or dict of str -> str models: Which models to plot. Can either be a list of models or a dict with mapping from model id to model name for display.
    :param bool snv: Whether to plot SNVs or Indels
    :param int cut: Bin cut on the entire data to highlight
    :param dict of str -> str colors: Optional colors to use (model name -> desired color)
    :return: Plots of the score distributions
    :rtype: Tabs
    """

    if not isinstance(models, dict):
        models = {m: m for m in models}

    if colors is None:
        colors = {m_name: "#033649" for m_name in models.values()}

    tabs = []
    for model_id, model_name in models.items():
        if model_id in ['vqsr', 'cnn', 'rf_2.0.2', 'rf_2.0.2_beta']:
            ht = hl.read_table(score_ranking_path(data_type, model_id, binned=False))
        else:
            ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=model_id))

        ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]), keep=snv)
        binned_ht = hl.read_table(score_ranking_path(data_type, model_id, binned=True))
        binned_ht = binned_ht.filter(binned_ht.snv, keep=snv)

        cut_value = binned_ht.aggregate(hl.agg.filter((binned_ht.bin == cut) & (binned_ht.rank_id == 'rank'), hl.agg.min(binned_ht.min_score)))

        min_score, max_score = (-20, 30) if model_id == 'vqsr' else (0.0, 1.0)
        agg_values = ht.aggregate(hl.struct(
            score_hist=[hl.agg.hist(ht.score, min_score, max_score, 100),
                        hl.agg.filter(ht.ac > 0, hl.agg.hist( ht.score, min_score, max_score, 100))],
            adj_counts=hl.agg.filter(ht.ac > 0, hl.agg.counter( ht.score >= cut_value))
        ))
        score_hist = agg_values.score_hist
        adj_cut = '{0:.2f}'.format(100 * agg_values.adj_counts[True] / (agg_values.adj_counts[True] + agg_values.adj_counts[False]))

        rows = []
        x_range = DataRange1d()
        y_range = [DataRange1d(), DataRange1d()]
        for adj in [False, True]:
            title = '{0}, {1} cut (score = {2:.2f})'.format('Adj' if adj else 'All', adj_cut if adj else cut, cut_value)
            p = plot_hail_hist(score_hist[adj], title=title + "\n", fill_color=colors[model_name])
            p.add_layout(Span(location=cut_value, dimension='height', line_color='red', line_dash='dashed'))
            p.x_range = x_range
            p.y_range = y_range[0]
            set_plots_defaults(p)

            p_cumul = plot_hail_hist_cumulative(score_hist[adj], title=title + ', cumulative', line_color=colors[model_name])
            p_cumul.add_layout(Span(location=cut_value, dimension='height', line_color='red', line_dash='dashed'))
            p_cumul.x_range = x_range
            p_cumul.y_range = y_range[1]
            set_plots_defaults(p_cumul)

            rows.append([p, p_cumul])

        tabs.append(Panel(child=gridplot(rows), title=model_name))
    return Tabs(tabs=tabs)


def get_binned_concordance_pd(data_type: str, truth_samples: List[str], models: Union[Dict[str, str], List[str]]) -> pd.DataFrame:
    """
    Creates a pandas DF containing the binned concordance results for all given truth samples / models.

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str truth_samples: List of truth samples to include
    :param list of str or dict of str -> str models: Models to include. Either a list of the model ids, or a dict with model id -> model name for display
    :return: Pandas dataframe with binned concordance results
    :rtype: DataFrame
    """

    def get_binned_concordance_ht(data_type: str, truth_samples: List[str], models: Dict[str, str]) -> hl.Table:
        """
        Combines binned concordance results for multiple truth samples and/or models into a single Table.
        """
        hts = []
        for truth_sample in truth_samples:
            for model_id, model_name in models.items():
                ht = hl.read_table(binned_concordance_path(data_type, truth_sample, model_id))
                ht = ht.annotate(truth_sample=truth_sample, model=model_name)
                hts.append(ht)

        return hts[0].union(*hts[1:])

    def compute_cumul_metrics(df: pd.DataFrame) -> pd.DataFrame:
        """
        Computes cumulative metrics on a pandas DF.
        """
        df = df.sort_values(by=['bin'])
        df['cum_tp'] = df['tp'].cumsum()
        df['cum_fp'] = df['fp'].cumsum()
        total_pos = df['tp'].sum() + df['fn'].sum()
        total_neg = df['fp'].sum()
        df['cum_tn'] = total_neg - df['cum_fp']
        df['cum_fn'] = total_pos - df['cum_tp']
        df['precision'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fp'])
        df['recall'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fn'])
        df['cum_alleles'] = df['n_alleles'].cumsum()
        return df[['bin', 'min_score', 'max_score', 'n_alleles', 'tp', 'fp', 'fn', 'cum_alleles', 'cum_tp', 'cum_fp', 'cum_fn', 'cum_tn', 'precision', 'recall']]

    if not isinstance(models, dict):
        models = {m: m for m in models}

    df = get_binned_concordance_ht(data_type, truth_samples, models).to_pandas()
    df = df.groupby(['rank_name', 'truth_sample', 'model', 'snv']).apply(compute_cumul_metrics)
    return df.fillna(-1).groupby(['rank_name', 'truth_sample', 'model', 'snv'])


def plot_concordance_pr(
        pr_df: pd.DataFrame,
        snv: bool,
        colors: Dict[str, str] = None,
        size_prop: str = None,
        bins_to_label: List[int] = None
) -> Column:
    """
    Generates plots showing Precision/Recall curves for truth samples:
    Two tabs:
    - One displaying the PR curve with ranking computed on the entire data
    - One displaying the PR curve with ranking computed on the  truth sample only

    Within each tab, a row of n_truth_samples.

    The input to this function should come out of the `get_binned_concordance_pd` function, which creates
    a DataFrame containing the necessary metris for PR plotting and is grouped by 'rank_name', 'truth_sample', 'model' and 'snv'.

    :param DataFrame pr_df: Input Dataframe
    :param bool snv: Whether to plot SNVs or Indels
    :param dict of str -> str colors: Optional colors to use (model name -> desired color)
    :param str size_prop: Either 'radius' or 'area' can be specified. If either is specified, the points will be sized proportionally to the amount of data in that point.
    :param list of int bins_to_label: Bins to label
    :return: Bokeh grid of plots
    :rtype: Tabs
    """

    if colors is None:
        # Get a palette automatically
        models = sorted(list(set([g[2] for g in pr_df.groups])))
        palette = d3['Category10'][max(3, len(models))]
        colors = {model: palette[i] for i, model in enumerate(models)}

    hover = HoverTool(tooltips=[
        ("model", "@model"),
        ("bin", "@bin"),
        ("score (min, max)", "(@min_score, @max_score)"),
        ('n_alleles', '@n_alleles'),
        ('cum_alleles', '@cum_alleles'),
        ("data (x,y)", "($x, $y)")
    ])

    tabs = []
    for rank in ['truth_sample_rank', 'global_rank']:

        plot_row = []
        for truth_sample in set([g[1] for g in pr_df.groups]):

            p = figure(title=truth_sample[0].upper() + truth_sample[1:],
                       x_axis_label='Recall',
                       y_axis_label='Precision',
                       tools=[hover] + [tool for tool in TOOLS.split(',') if tool != 'hover'])
            p.xaxis[0].formatter = NumeralTickFormatter(format="0%")
            p.yaxis[0].formatter = NumeralTickFormatter(format="0.0%")

            circles = []
            for model in set([g[2] for g in pr_df.groups]):
                data = pr_df.get_group((rank, truth_sample, model, snv)).copy()
                data['model'] = [model] * len(data)
                data['size'] = get_point_size_col(data['n_alleles'], size_prop)
                source = ColumnDataSource(data)
                circles.append((model, [p.circle('recall',
                         'precision',
                         size='size',
                         color=colors[model], source=source)]))
                if bins_to_label is not None:
                    label_data = data.loc[data.bin.isin(bins_to_label)].copy()
                    label_data['x_offset'] = label_data['recall'] + 0.025
                    label_data['y_offset'] = label_data['precision']
                    label_data['bin_str'] = [str(int(t)) for t in label_data['bin']]
                    label_source = ColumnDataSource(label_data)
                    p.add_layout(
                        LabelSet(x='x_offset',
                                 y='precision',
                                 text='bin_str',
                                 text_color=colors[model],
                                 source=label_source)
                    )
                    p.multi_line(
                        xs=[[x, x + 0.05] for x in label_data.recall],
                        ys=[[y, y] for y in label_data.precision],
                        color=colors[model]
                    )

            legend = Legend(items=circles, orientation='horizontal', location=(0, 0), click_policy="hide")
            p.add_layout(legend, 'above')
            set_plots_defaults(p)
            plot_row.append(p)

        tabs.append(Panel(child=Row(children=plot_row), title=rank))

    return Tabs(tabs=tabs)

