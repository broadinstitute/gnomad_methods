from gnomad_hail import *
from gnomad_hail.resources.variant_qc import *
from gnomad_hail.utils.plotting import *


def get_binned_models_pd(data_type: str, models: List[str], contigs: Set[str] = None) -> pd.DataFrame:
    """
    Creates a single DataFrame with all desired models binned and ready for plotting.

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str models: Models to load
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
            max_score=hl.agg.min(ht.max_score),
            **{x: hl.agg.sum(ht[x]) for x in ht.row_value if x not in ['min_score', 'max_score']}
        )

    hts = [
        aggregate_contig(hl.read_table(score_ranking_path(data_type, model, binned=True)), contigs)
            .annotate(model=model) for model in models
    ]

    ht = hts.pop()
    if hts:
        ht = ht.union(*hts)
    return ht.to_pandas()


def plot_metric(df: pd.DataFrame,
                y_name: str,
                cols: List[str],
                y_fun: Callable[[pd.Series], Union[float, int]] = lambda x: x,
                plot_all: bool = True,
                plot_bi_allelics: bool = True,
                plot_singletons: bool = True,
                plot_bi_allelic_singletons: bool = True,
                plot_release: bool = False,
                colors: Dict[str, str] = None,
                link_cumul_y: bool = True,
                legend_position: str = 'top_right') -> Tabs:
    """
    Generic function for generating QC metric plots using a plotting-ready DataFrame (obtained from `get_binned_models_pd`)
    DataFrame needs to have a `rank_id` column, a `bin` column and a `model` column (contains the model name and needs to be added to binned table(s))

    This function generates scatter plots with the metric bin on x-axis and a user-defined function on the y-axis.
    The data for the y-axis function needs to from the columns specified in `cols`. The function is specified with the `y_fun` argument and data columns are access as a list.
    As an example, plotting Transition to transversion ratio is done as follows:
    ```
    plot_metric(snvs, 'Ti/Tv', ['n_ti', 'n_tv'], y_fun=lambda x: x[0]/x[1], colors=colors)

    ```
    In this command, `x[0]` correspond to the  first column selected (`'n_ti'`)  and `x[1]` to the second (`'n_tv'`).


    This function plots a tab for each of the plot condition(s) selected: all, bi-allelics, bi-allelic singletons.
    Within each tab, each row contains a non-cumulative and a cumulative plot of the bins / values.
    If `plot_release` is set, then an extra row is added plotting release-sites only (variants in release samples where AC_ADJ>0). The bin for these sites is computed based on relase variants only.

    :param pd.DataFrame df: Input data
    :param str y_name: Name of the metric plotted on the y-axis
    :param list of str cols: Columns used to compute the metric plotted
    :param callable y_fun: Function to apply to the columns to generate the metric
    :param bool plot_all: Whether to plot a tab with all variants
    :param bool plot_bi_allelics: Whether to plot a tab with bi-allelic variants only
    :param bool plot_singletons: Whether to plot a tab with singleton variants only
    :param bool plot_bi_allelic_singletons:  Whether to plot a tab with bi-allelic singleton variants only
    :param bool plot_release: Whether to plot additional rows with release-only variants
    :param dict of str -> str colors: Mapping of model name -> color
    :param bool link_cumul_y: If set, y-axes of cumulative and non-cumulative plots are linked
    :param str legend_position: Legend position in the plot
    :return: Plot
    :rtype: Tabs
    """

    def get_row(df: pd.DataFrame, y_name: str, cols: List[str], y_fun: Callable[[pd.Series], Union[float, int]], title: str, link_cumul_y: bool, legend_position: str) -> Row:
        """
        Generates a single row
        """

        def get_plot(df: pd.DataFrame, y_name: str, y_col_name: str, title: str, data_ranges: Tuple[DataRange1d, DataRange1d], legend_position: str) -> Plot:
            """
            Generates a single plot panel
            """
            data_source = ColumnDataSource(df)
            p = figure(
                title=title,
                x_axis_label='bin',
                y_axis_label=y_name,
                tools="save,pan,box_zoom,reset,wheel_zoom,box_select,lasso_select,help,hover")
            p.x_range = data_ranges[0]
            p.y_range = data_ranges[1]
            p.circle('bin', y_col_name, legend='model', color='color', source=data_source)
            p.select_one(HoverTool).tooltips = [('model', '@model'),
                                                ('bin', '@bin'),
                                                (y_name, f'@{y_col_name}'),
                                                ('min_score', '@min_score'),
                                                ('max_score', '@max_score')
                                                ] + [(col, f'@{col}') for col in cols]
            p.legend.location = legend_position
            return p

        df['non_cumul'] = df[cols].apply(y_fun, axis=1)
        for col in cols:
            df[f'{col}_cumul'] = df.groupby('model').aggregate(np.cumsum)[col]

        df['cumul'] = df[[f'{col}_cumul' for col in cols]].apply(y_fun, axis=1)
        non_cumul_data_ranges = (DataRange1d(), DataRange1d())
        cumul_data_ranges = non_cumul_data_ranges if link_cumul_y else (non_cumul_data_ranges[0], DataRange1d())

        return Row(get_plot(df, y_name, 'non_cumul', title, non_cumul_data_ranges, legend_position),
                   get_plot(df, y_name, 'cumul', title + ', cumulative', cumul_data_ranges, legend_position))

    def prepare_pd(df: pd.DataFrame, cols: List[str], colors: Dict[str, str] = {}):
        """
        Groups a pandas DataFrame by model and bin while keeping relevant columns only. Adds a color column.
        """
        df = df.groupby(['model', 'bin']).agg({**{col: np.sum for col in cols},
                                               'min_score': np.min, 'max_score': np.max})
        df = df.reset_index()
        df['color'] = [colors.get(x, 'gray') for x in df['model']]
        return df

    colors = colors if colors is not None else {}
    tabs = []
    release_strats = ['', 'release_'] if plot_release else ['']

    if plot_all:
        children = []
        for release in release_strats:
            title = '{0}, {1}'.format(y_name, 'release' if release else 'all')
            children.append(get_row(prepare_pd(df.loc[df.rank_id == f'{release}rank'], cols, colors),
                                    y_name, cols, y_fun, title, link_cumul_y, legend_position))

        tabs.append(Panel(child=Column(children=children), title='All'))

    if plot_bi_allelics:
        children = []
        for release in release_strats:
            for biallelic_rank in ['', 'biallelic_']:
                title = '{0}, bi-allelic ({1} rank)'.format(y_name, 'overall' if not release and not biallelic_rank else f'{release[:-1]} {biallelic_rank[:-1]}')
                children.append(get_row(prepare_pd(df.loc[df.bi_allelic & (df.rank_id == f'{release}{biallelic_rank}rank')], cols, colors),
                                        y_name, cols, y_fun, title, link_cumul_y, legend_position))

        tabs.append(Panel(child=Column(children=children), title='Bi-allelic'))

    if plot_singletons:
        children = []
        for release in release_strats:
            for singleton_rank in ['', 'singleton_']:
                title = '{0}, singletons ({1} rank)'.format(y_name, 'overall' if not release and not singleton_rank else f'{release[:-1]} {singleton_rank[:-1]}')
                children.append(get_row(prepare_pd(df.loc[df.singleton & (df.rank_id == f'{release}{singleton_rank}rank')], cols, colors),
                                        y_name, cols, y_fun, title, link_cumul_y, legend_position))

        tabs.append(Panel(child=Column(children=children), title='Singletons'))

    if plot_bi_allelic_singletons:
        children = []
        for release in release_strats:
            for bisingleton_rank in ['', 'biallelic_singleton_']:
                title = '{0}, singletons ({1} rank)'.format(y_name, 'overall' if not release and not bisingleton_rank else f'{release[:-1]} {bisingleton_rank[:-1].replace("_", " ")}')
                children.append(get_row(prepare_pd(df.loc[df.bi_allelic & df.singleton & (df.rank_id == f'{release}{bisingleton_rank}rank')], cols, colors),
                                        y_name, cols, y_fun, title, link_cumul_y, legend_position))

        tabs.append(Panel(child=Column(children=children), title='Bi-allelic singletons'))

    return Tabs(tabs=tabs)


def plot_score_distributions(data_type, models: List[str], snv: bool, cut: int) -> Tabs:
    """
    Generates plots of model scores distributions:
    One tab per model.
    Within each tab, there is 2x2 grid of plots:
    - One row showing the score distribution across the entire data
    - One row showing the score distribution across the release data only (release_sample_AC_ADJ > 0)
    - One column showing the histogram of the score
    - One column showing the normalized cumulative histogram of the score

    Cutoff is highlighted by a dashed red line

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str models: Which models to plot
    :param bool snv: Whether to plot SNVs or Indels
    :param int cut: Bin cut on the entire data to highlight
    :return: Plots of the score distributions
    :rtype: Tabs
    """

    tabs = []
    for model in models:
        if model in ['vqsr', 'cnn', 'rf_2.0.2', 'rf_2.0.2_beta']:
            ht = hl.read_table(score_ranking_path(data_type, model, binned=False))
        else:
            ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=model))

        ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]), keep=snv)
        binned_ht = hl.read_table(score_ranking_path(data_type, model, binned=True))
        binned_ht = binned_ht.filter(binned_ht.snv, keep=snv)

        cut_value = binned_ht.aggregate(hl.agg.min(hl.agg.filter((binned_ht.bin == cut) & (binned_ht.rank_id == 'rank'), binned_ht.min_score)))

        min_score, max_score = (-20, 30) if model == 'vqsr' else (0.0, 1.0)
        agg_values = ht.aggregate(hl.struct(
            score_hist=[hl.agg.hist(ht.score, min_score, max_score, 100),
                        hl.agg.hist(hl.agg.filter(ht.ac > 0, ht.score), min_score, max_score, 100)],
            release_counts=hl.agg.counter(hl.agg.filter(ht.ac > 0, ht.score >= cut_value))
        ))
        score_hist = agg_values.score_hist
        release_cut = '{0:.2f}'.format(100 * agg_values.release_counts[True] / (agg_values.release_counts[True] + agg_values.release_counts[False]))

        rows = []
        x_range = DataRange1d()
        y_range = [DataRange1d(), DataRange1d()]
        for release in [False, True]:
            title = '{0}, {1} cut (score = {2:.2f})'.format('Release' if release else 'All', release_cut if release else cut, cut_value)
            p = plot_hail_hist(score_hist[release], title=title)
            p.line(x=[cut_value, cut_value], y=[0, max(score_hist[release].bin_freq)], color='red', line_dash='dashed')
            p.x_range = x_range
            p.y_range = y_range[0]
            p_cumul = plot_hail_hist_cumulative(score_hist[release], title=title + ' cumulative')
            p_cumul.line(x=[cut_value, cut_value], y=[0.0, 1.0], color='red', line_dash='dashed')
            p_cumul.x_range = x_range
            p_cumul.y_range = y_range[1]

            rows.append([p, p_cumul])

        tabs.append(Panel(child=gridplot(rows), title=model))
    return Tabs(tabs=tabs)


def get_binned_concordance_pd(data_type: str, truth_samples: List[str], models: List[str]) -> pd.DataFrame:
    """
    Creates a pandas DF containing the binned concordance results for all given truth samples / models.

    :param str data_type: One of 'exomes' or 'genomes'
    :param list of str truth_samples: List of truth samples to include
    :param list of str models: List of models to include
    :return: Pandas dataframe with binned concordance results
    :rtype: DataFrame
    """

    def get_binned_concordance_ht(data_type: str, truth_samples: List[str], models: List[str]) -> hl.Table:
        """
        Combines binned concordance results for multiple truth samples and/or models into a single Table.
        """
        hts = []
        for truth_sample in truth_samples:
            for metric in models:
                ht = hl.read_table(binned_concordance_path(data_type, truth_sample, metric))
                ht = ht.annotate(truth_sample=truth_sample, metric=metric)
                hts.append(ht)

        return hts[0].union(*hts[1:])

    def compute_cumul_metrics(df: pd.DataFrame) -> pd.DataFrame:
        """
        Computes cumulative metrics on a pandas DF.
        """
        df = df.sort_values(by=['bin'])
        df['cum_tp'] = df['tp'].cumsum()
        df['cum_fp'] = df['fp'].cumsum()
        total_tps = df['tp'].sum() + df['fn'].sum()
        total_fps = df['fp'].sum()
        df['cum_tn'] = total_fps - df['cum_fp']
        df['cum_fn'] = total_tps - df['cum_tp']
        df['precision'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fp'])
        df['recall'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fn'])
        df['cum_alleles'] = df['n_alleles'].cumsum()
        return df[['bin', 'min_score', 'max_score', 'n_alleles', 'tp', 'fp', 'fn', 'cum_alleles', 'cum_tp', 'cum_fp', 'cum_fn', 'cum_tn', 'precision', 'recall']]

    df = get_binned_concordance_ht(data_type, truth_samples, models).to_pandas()
    df = df.groupby(['rank_name', 'truth_sample', 'metric', 'snv']).apply(compute_cumul_metrics)
    return df.fillna(-1).groupby(['rank_name', 'truth_sample', 'metric', 'snv'])


def plot_concordance_pr(
        pr_df: pd.DataFrame,
        data_type: str,
        plot_colors: Dict[str, str] = None,
        legend_text: Dict[str, str] = None,
        size_prop: str = None,
        bins_to_label: Dict[str, List[int]] = None
) -> Column:
    """
    Generates plots showing Precision/Recall curves for truth samples:
    Two tabs:
    - One displaying the PR curve with ranking computed on the entire data
    - One displaying the PR curve with ranking computed on the  truth sample only

    Within each tab, a grid of n_truth_samples x 2 plots:
    - One row per truth sample
    - One columns per allele type (SNVs, Indels)

    :param DataFrame pr_df: Input Dataframe
    :param str data_type: One of 'exomes' or 'genomes' -- use in the title only
    :param dict of str -> str plot_colors: Optional colors to use (model name -> desired color)
    :param dict of str -> str legend_text: Optional mapping of model to model name
    :param str size_prop: Either 'size' or 'area' can be specified. If either is specified, the points will be sized proportionally to the amount of data in that point.
    :param dict of str -> list of int bins_to_label: Bins to label for each variant type ('snv', 'indel')
    :return Bokeh grid of plots
    :rtype Tabs
    """

    if plot_colors is None:
        # Get a palette automatically
        from bokeh.palettes import d3
        models = sorted(list(set([g[2] for g in pr_df.groups])))
        palette = d3['Category10'][max(3, len(models))]
        plot_colors = {model: palette[i] for i, model in enumerate(models)}

    if legend_text is None:
        legend_text = {}

    hover = HoverTool(tooltips=[
        ("bin", "@bin"),
        ("score (min, max)", "(@min_score, @max_score)"),
        ('n_alleles', '@n_alleles'),
        ('cum_alleles', '@cum_alleles'),
        ("data (x,y)", "($x, $y)"),

    ])

    tabs = []
    for rank in ['global_rank', 'truth_sample_rank']:

        plot_grid = []
        for truth_sample in set([g[1] for g in pr_df.groups]):
            plot_rows = []
            for snv in [True, False]:
                variant_type = 'snv' if snv else 'indel'
                p = figure(title='{} {}, {}'.format(data_type, truth_sample, variant_type),
                           x_axis_label='Recall',
                           y_axis_label='Precision',
                           tools=[hover] + [tool for tool in TOOLS.split(',') if tool != 'hover'])
                p.xaxis[0].formatter = NumeralTickFormatter(format="0%")
                p.yaxis[0].formatter = NumeralTickFormatter(format="0.0%")
                p.title.text_font_size = "20pt"
                p.axis.axis_label_text_font_size = "20pt"
                p.axis.axis_label_text_font_style = "normal"
                p.axis.major_label_text_font_size = "14pt"
                for model in set([g[2] for g in pr_df.groups]):
                    data = pr_df.get_group((rank, truth_sample, model, snv))
                    total_alleles = np.sum(data['n_alleles'])
                    if size_prop is None:
                        data['radius'] = 4
                    elif size_prop == 'radius':
                        data['radius'] = 4 * (data['n_alleles'] / total_alleles)
                    elif size_prop == 'area':
                        data['radius'] = 8 * np.sqrt((data['n_alleles'] / total_alleles) / np.pi)
                    else:
                        raise ValueError(f"{size_prop} is not a supported value for argument `size_prop`")
                    source = ColumnDataSource(data)
                    p.circle('recall',
                             'precision',
                             radius='radius',
                             radius_units='screen',
                             legend=legend_text.get(model, model),
                             color=plot_colors[model], source=source)
                    if bins_to_label is not None:
                        label_data = data.loc[data.bin.isin(bins_to_label[variant_type])]
                        label_data['x_offset'] = label_data['recall'] + 0.025
                        label_data['y_offset'] = label_data['precision']
                        label_data['bin_str'] = [str(int(t)) for t in label_data['bin']]
                        label_source = ColumnDataSource(label_data)
                        p.add_layout(
                            LabelSet(x='x_offset',
                                     y='precision',
                                     text='bin_str',
                                     text_color=plot_colors[model],
                                     source=label_source)
                        )
                        p.multi_line(
                            xs=[[x, x + 0.05] for x in label_data.recall],
                            ys=[[y, y] for y in label_data.precision],
                            color=plot_colors[model]
                        )

                p.legend.location = 'bottom_left'
                p.legend.label_text_font_size = "14pt"
                plot_rows.append(p)

            plot_grid.append(plot_rows)

        tabs.append(Panel(child=gridplot(plot_grid), title=rank))

    return Tabs(tabs=tabs)
