
import hail as hl
import numpy as np
from ipywidgets import interact
import math
import pandas as pd
from collections import OrderedDict
import json
from .constants import *

import bokeh
from bokeh.layouts import gridplot, row, widgetbox
from bokeh.plotting import figure, show, output_file
from bokeh.io import output_notebook, push_notebook, export_png
from bokeh.models.widgets import Tabs, Panel
from bokeh.palettes import *
from bokeh.models import *
from typing import *
from bokeh.plotting.helpers import stack
from bokeh.transform import factor_cmap

# Setting some defaults for Table.show
if 'old_show' not in dir():
    old_show = hl.Table.show

    def new_show(t, n=10, width=170, truncate=40, types=True):
        old_show(t, n, width, truncate, types)
    hl.Table.show = new_show

TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

COLOR_CPG = '#2E9FFE'
COLOR_TI = '#458B00'
COLOR_TV = '#EA4444'
ACG = '#422efe'
CCG = '#2e73fe'
GCG = '#2ec3fe'
TCG = '#2efef7'
COLOR_METHYLATED_CPG = '#2E37FE'

variant_type_colors = {
    'CpG': COLOR_CPG,
    'non-CpG': COLOR_TI,
    'non-CpG transition': COLOR_TI,
    'transversion': COLOR_TV,
    'ACG': ACG, 'CCG': CCG, 'GCG': GCG, 'TCG': TCG,
    'methylated CpG': COLOR_METHYLATED_CPG
}

COLOR_SYN = '#AAAAAA'
COLOR_MIS = '#FF6103'
COLOR_LOF = '#9D1309'

variant_annotation_colors = {x: COLOR_LOF for x in CSQ_CODING_HIGH_IMPACT}
variant_annotation_colors.update({x: COLOR_MIS for x in CSQ_CODING_MEDIUM_IMPACT})
variant_annotation_colors.update({x: COLOR_SYN for x in CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING})
variant_annotation_colors.update({
    'stop_lost': COLOR_MIS,
    'splice_region_variant': COLOR_SYN,
    'start_lost': COLOR_SYN,
    'Synonymous': COLOR_SYN,
    'Missense': COLOR_MIS,
    'LoF': COLOR_LOF,
    'PTV': COLOR_LOF
})

COLOR_EXAC = '#4682B4'
COLOR_GNOMAD = '#73AB3D'

dataset_colors = {
    'ExAC': COLOR_EXAC,
    'gnomAD': COLOR_GNOMAD
}

COLOR_AMR = '#ED1E24'
COLOR_EUR = '#6AA5CD'
COLOR_AFR = '#941494'
COLOR_SAS = '#FF9912'
COLOR_EAS = '#108C44'
COLOR_OTH = '#ABB9B9'
COLOR_MDE = '#33CC33'
COLOR_ASJ = '#FF7F50'
COLOR_NFE = COLOR_EUR
COLOR_FIN = '#002F6C'


pop_colors = {'afr': COLOR_AFR,
              'amr': COLOR_AMR,
              'eas': COLOR_EAS,
              'fin': COLOR_FIN,
              'eur': COLOR_NFE,
              'nfe': COLOR_NFE,
              'oth': COLOR_OTH,
              'sas': COLOR_SAS,
              'mde': COLOR_MDE,
              'asj': COLOR_ASJ,
              'uniform': 'pink',
              'consanguineous': 'pink',
              'sas_non_consang': 'orange',
              'exac': 'gray',
              'est': 'black',
              'bg': '#66C2A5',
              'de': 'black',
              'ee': '#4891D9',
              'es': '#FFC400',
              'fi': COLOR_FIN,
              'neu': '#C60C30',
              'seu': '#009246',
              'gb': '#C60C30',
              'it': '#009246',
              'se': 'purple',
              'cn': '#FFC400',
              'kr': '#4891D9',
              'nka': '#009246',
              'hk': '#C60C30',
              'sg': 'darkred',
              'tw': '#009246',
              'jp': '#BC002D',
              'unk': COLOR_OTH,
              '': COLOR_OTH}
pop_names = {'oth': 'Other',
             'afr': 'African',
             'amr': 'Latino',
             'eas': 'East Asian',
             'fin': 'Finnish',
             'eur': 'European',
             'nfe': 'European',
             'sas': 'South Asian',
             'mde': 'Middle Eastern',
             'asj': 'Ashkenazi Jewish',
             'uniform': 'Uniform',
             'sas_non_consang': 'South Asian (F < 0.05)',
             'consanguineous': 'South Asian (F > 0.05)',
             'exac': 'ExAC',
             'bg': 'Bulgaria (Eastern Europe)',
             'de': 'Germany',
             'ee': 'Estonia',
             'es': 'Spain',
             'fi': 'Finland',
             'gb': 'UK',
             'nwe': 'North-Western Europe',
             'seu': 'Southern Europe',
             'it': 'Italy',
             'se': 'Sweden',
             'cn': 'China',
             'kr': 'Korea',
             'hk': 'Hong Kong',
             'sg': 'Singapore',
             'tw': 'Taiwan',
             'jp': 'Japan',
             'oea': 'Other Asian',
             'oeu': 'Other European',
             't2d': 'T2DGenes',
             'jhs': 'JHS',
             'biome': 'BioMe',
             'unk': 'Unknown'}


def plot_hail_hist(hist_data: hl.Struct,
                   title: str = 'Plot',
                   log: bool = False,
                   fill_color: str = "#033649",
                   outlier_fill_color: str = "#036564",
                   line_color: str = '#033649',
                   hover_mode: str = 'vline') -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))

    :param Struct hist_data: Data to plot
    :param str title: Plot title
    :param bool log: Whether the y-axis should be log
    :param str fill_color: Color to fill the histogram bars that fall within the hist boundaries
    :param outlier_fill_color: Color to fill the histogram bars that fall outside the hist boundaries
    :param str line_color: Color of the lines around the histogram bars
    :param str hover_mode: Hover mode; one of 'vline' (default), 'mouse' or 'hline'
    :return: Histogram plot
    :rtype: Figure
    """
    low = int(log)
    distance = abs(hist_data.bin_edges[0] - hist_data.bin_edges[1])
    num_data_points = sum(hist_data.bin_freq)
    p = figure(title=title, y_axis_type="log", tools=TOOLS) if log else figure(title=title, tools=TOOLS)
    p.add_layout(Title(text=f'({num_data_points:,} data points)'), 'above')

    top = [x + low for x in hist_data.bin_freq]
    left = hist_data.bin_edges[:-1]
    right = hist_data.bin_edges[1:]
    colors = [fill_color] * len(hist_data.bin_freq)
    if hist_data.n_smaller > 0:
        top.insert(0, hist_data.n_smaller + low)
        left.insert(0, hist_data.bin_edges[0] - distance)
        right.insert(0, hist_data.bin_edges[0])
        colors.insert(0, outlier_fill_color)
    if hist_data.n_larger > 0:
        top.append(hist_data.n_larger + low)
        left.append(hist_data.bin_edges[-1])
        right.append(hist_data.bin_edges[-1] + distance)
        colors.append(outlier_fill_color)

    hist_source = ColumnDataSource({'top': top,
                                    'bottom': [low] * len(top),
                                    'left': left,
                                    'right': right,
                                    'color': colors
                                    })

    p.select_one(HoverTool).tooltips = [("bin", "$index"), ("bin_edges", "(@left, @right)"), ("freq", "@top")]
    p.select_one(HoverTool).mode = hover_mode
    hist_data.bin_freq = [x + low for x in hist_data.bin_freq]
    p.quad(top='top', bottom='bottom',
           left='left', right='right', source=hist_source, fill_color='color', line_color=line_color)

    return p


def plot_hail_hist_cumulative(hist_data: hl.Struct, title: str = 'Plot', normalize: bool = True,
                              line_color: str = "#036564", line_width: int = 3, log: bool = False, hover_mode: str = 'mouse') -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))

    :param Struct hist_data: Data to plot
    :param str title: Plot title
    :param bool normalize: Whether to normalize the data (0,1)
    :param str line_color: Color of the line
    :param int line_width: Width of the line
    :param bool log: Whether the y-axis should be log
    :param str hover_mode: Hover mode; one of 'mouse' (default), 'vline' or 'hline'
    :return: Histogram plot
    :rtype: Figure
    """
    cumulative_data = np.cumsum(hist_data.bin_freq) + hist_data.n_smaller
    np.append(cumulative_data, [cumulative_data[-1] + hist_data.n_larger])
    num_data_points = max(cumulative_data)

    if normalize: cumulative_data = cumulative_data / num_data_points
    p = figure(title=title, y_axis_type="log", tools=TOOLS) if log else figure(title=title, tools=TOOLS)
    p.add_layout(Title(text=f'({num_data_points:,} data points)'), 'above')
    p.select_one(HoverTool).tooltips = [("index", "$index"), ("(x,y)", "(@x, @y)")]
    p.select_one(HoverTool).mode = hover_mode
    data_source = ColumnDataSource({'x': hist_data.bin_edges[:-1], 'y': cumulative_data})
    p.line(x='x', y='y', line_color=line_color, line_width=line_width, source=data_source)
    return p


def plot_hail_hist_both(hist_data: hl.Struct, title: str, normalize: bool = True, log: bool = False):
    p1 = plot_hail_hist(hist_data, title, log)
    p2 = plot_hail_hist_cumulative(hist_data, f'{title} (Cumulative)', normalize, log=log)
    return Tabs(tabs=[Panel(child=p1, title='Raw'), Panel(child=p2, title='Cumulative')])


def set_font_size(p, font_size: str = "12pt"):
    p.title.text_font_size = font_size
    p.legend.label_text_font_size = font_size
    p.xaxis.axis_label_text_font_size = font_size
    p.yaxis.axis_label_text_font_size = font_size
    p.xaxis.major_label_text_font_size = font_size
    p.yaxis.major_label_text_font_size = font_size
    if hasattr(p.xaxis, 'group_text_font_size'):
        p.xaxis.group_text_font_size = font_size
    return p


def linear_and_log_tabs(plot_func: Callable, **kwargs) -> Tabs:
    panels = []
    for axis_type in ["linear", "log"]:
        fig = plot_func(**kwargs, axis_type=axis_type)
        panel = Panel(child=fig, title=axis_type)
        panels.append(panel)

    return Tabs(tabs=panels)


def plot_hail_file_metadata(t_path: str) -> Optional[Union[Grid, Tabs, bokeh.plotting.Figure]]:
    """
    Takes path to hail Table or MatrixTable (gs://bucket/path/hail.mt), outputs Grid or Tabs, respectively
    Or if an unordered Table is provided, a Figure with file sizes is output
    If metadata file or rows directory is missing, returns None
    """
    panel_size = 600
    subpanel_size = 150

    files = hl.hadoop_ls(t_path)
    rows_file = [x['path'] for x in files if x['path'].endswith('rows')]
    entries_file = [x['path'] for x in files if x['path'].endswith('entries')]
    # cols_file = [x['path'] for x in files if x['path'].endswith('cols')]
    success_file = [x['modification_time'] for x in files if x['path'].endswith('SUCCESS')]

    data_type = 'Table'

    metadata_file = [x['path'] for x in files if x['path'].endswith('metadata.json.gz')]
    if not metadata_file:
        warnings.warn('No metadata file found. Exiting...')
        return None

    with hl.hadoop_open(metadata_file[0], 'rb') as f:
        overall_meta = json.loads(f.read())
        rows_per_partition = overall_meta['components']['partition_counts']['counts']

    if not rows_file:
        warnings.warn('No rows directory found. Exiting...')
        return None
    rows_files = hl.hadoop_ls(rows_file[0])

    if entries_file:
        data_type = 'MatrixTable'
        rows_file = [x['path'] for x in rows_files if x['path'].endswith('rows')]
        rows_files = hl.hadoop_ls(rows_file[0])
    row_partition_bounds, row_file_sizes = get_rows_data(rows_files)

    total_file_size, row_file_sizes, row_scale = scale_file_sizes(row_file_sizes)

    if not row_partition_bounds:
        warnings.warn('Table is not partitioned. Only plotting file sizes')
        row_file_sizes_hist, row_file_sizes_edges = np.histogram(row_file_sizes, bins=50)
        p_file_size = figure(plot_width=panel_size, plot_height=panel_size)
        p_file_size.quad(right=row_file_sizes_hist, left=0, bottom=row_file_sizes_edges[:-1],
                         top=row_file_sizes_edges[1:], fill_color="#036564", line_color="#033649")
        p_file_size.yaxis.axis_label = f'File size ({row_scale}B)'
        return p_file_size

    all_data = {
        'partition_widths': [-1 if x[0] != x[2] else x[3] - x[1] for x in row_partition_bounds],
        'partition_bounds': [f'{x[0]}:{x[1]}-{x[2]}:{x[3]}' for x in row_partition_bounds],
        'spans_chromosome': ['Spans chromosomes' if x[0] != x[2] else 'Within chromosome' for x in row_partition_bounds],
        'row_file_sizes': row_file_sizes,
        'row_file_sizes_human': [f'{x:.1f} {row_scale}B' for x in row_file_sizes],
        'rows_per_partition': rows_per_partition,
        'index': list(range(len(rows_per_partition)))
    }

    if entries_file:
        entries_rows_files = hl.hadoop_ls(entries_file[0])
        entries_rows_file = [x['path'] for x in entries_rows_files if x['path'].endswith('rows')]
        if entries_rows_file:
            entries_files = hl.hadoop_ls(entries_rows_file[0])
            entry_partition_bounds, entry_file_sizes = get_rows_data(entries_files)
            total_entry_file_size, entry_file_sizes, entry_scale = scale_file_sizes(entry_file_sizes)
            all_data['entry_file_sizes'] = entry_file_sizes
            all_data['entry_file_sizes_human'] = [f'{x:.1f} {entry_scale}B' for x in row_file_sizes]

    title = f'{data_type}: {t_path}'

    msg = f"Rows: {sum(all_data['rows_per_partition']):,}<br/>Partitions: {len(all_data['rows_per_partition']):,}<br/>Size: {total_file_size}<br/>"
    if success_file[0]:
        msg += success_file[0]

    source = ColumnDataSource(pd.DataFrame(all_data))
    p = figure(tools=TOOLS, plot_width=panel_size, plot_height=panel_size)
    p.title.text = title
    p.xaxis.axis_label = 'Number of rows'
    p.yaxis.axis_label = f'File size ({row_scale}B)'
    color_map = factor_cmap('spans_chromosome', palette=Spectral8,
                            factors=list(set(all_data['spans_chromosome'])))
    p.scatter('rows_per_partition', 'row_file_sizes', color=color_map, legend='spans_chromosome', source=source)
    p.legend.location = 'bottom_right'
    p.select_one(HoverTool).tooltips = [(x, f'@{x}') for x in
                                        ('rows_per_partition', 'row_file_sizes_human', 'partition_bounds', 'index')]

    p_stats = Div(text=msg)
    p_rows_per_partition = figure(x_range=p.x_range, plot_width=panel_size, plot_height=subpanel_size)
    p_file_size = figure(y_range=p.y_range, plot_width=subpanel_size, plot_height=panel_size)

    rows_per_partition_hist, rows_per_partition_edges = np.histogram(all_data['rows_per_partition'], bins=50)
    p_rows_per_partition.quad(top=rows_per_partition_hist, bottom=0, left=rows_per_partition_edges[:-1],
                              right=rows_per_partition_edges[1:],
                              fill_color="#036564", line_color="#033649")
    row_file_sizes_hist, row_file_sizes_edges = np.histogram(all_data['row_file_sizes'], bins=50)
    p_file_size.quad(right=row_file_sizes_hist, left=0, bottom=row_file_sizes_edges[:-1],
                     top=row_file_sizes_edges[1:], fill_color="#036564", line_color="#033649")

    rows_grid = gridplot([[p_rows_per_partition, p_stats], [p, p_file_size]])

    if 'entry_file_sizes' in all_data:
        title = f'Statistics for {data_type}: {t_path}'

        msg = f"Rows: {sum(all_data['rows_per_partition']):,}<br/>Partitions: {len(all_data['rows_per_partition']):,}<br/>Size: {total_entry_file_size}<br/>"
        if success_file[0]:
            msg += success_file[0]

        source = ColumnDataSource(pd.DataFrame(all_data))
        panel_size = 600
        subpanel_size = 150
        p = figure(tools=TOOLS, plot_width=panel_size, plot_height=panel_size)
        p.title.text = title
        p.xaxis.axis_label = 'Number of rows'
        p.yaxis.axis_label = f'File size ({entry_scale}B)'
        color_map = factor_cmap('spans_chromosome', palette=Spectral8, factors=list(set(all_data['spans_chromosome'])))
        p.scatter('rows_per_partition', 'entry_file_sizes', color=color_map, legend='spans_chromosome', source=source)
        p.legend.location = 'bottom_right'
        p.select_one(HoverTool).tooltips = [(x, f'@{x}') for x in ('rows_per_partition', 'entry_file_sizes_human', 'partition_bounds', 'index')]

        p_stats = Div(text=msg)
        p_rows_per_partition = figure(x_range=p.x_range, plot_width=panel_size, plot_height=subpanel_size)
        p_rows_per_partition.quad(top=rows_per_partition_hist, bottom=0, left=rows_per_partition_edges[:-1],
                                  right=rows_per_partition_edges[1:],
                                  fill_color="#036564", line_color="#033649")
        p_file_size = figure(y_range=p.y_range, plot_width=subpanel_size, plot_height=panel_size)

        row_file_sizes_hist, row_file_sizes_edges = np.histogram(all_data['entry_file_sizes'], bins=50)
        p_file_size.quad(right=row_file_sizes_hist, left=0, bottom=row_file_sizes_edges[:-1],
                         top=row_file_sizes_edges[1:], fill_color="#036564", line_color="#033649")
        entries_grid = gridplot([[p_rows_per_partition, p_stats], [p, p_file_size]])

        return Tabs(tabs=[Panel(child=entries_grid, title='Entries'), Panel(child=rows_grid, title='Rows')])
    else:
        return rows_grid


def scale_file_sizes(file_sizes):
    min_file_size = min(file_sizes) * 1.1
    total_file_size = sum(file_sizes)
    all_scales = [
        ('T', 1e12),
        ('G', 1e9),
        ('M', 1e6),
        ('K', 1e3),
        ('', 1e0)
    ]
    for overall_scale, overall_factor in all_scales:
        if total_file_size > overall_factor:
            total_file_size /= overall_factor
            break
    for scale, factor in all_scales:
        if min_file_size > factor:
            file_sizes = [x / factor for x in file_sizes]
            break
    total_file_size = f'{total_file_size:.1f} {overall_scale}B'
    return total_file_size, file_sizes, scale


def get_rows_data(rows_files):
    file_sizes = []
    partition_bounds = []
    parts_file = [x['path'] for x in rows_files if x['path'].endswith('parts')]
    if parts_file:
        parts = hl.hadoop_ls(parts_file[0])
        for i, x in enumerate(parts):
            index = x['path'].split(f'{parts_file[0]}/part-')[1].split('-')[0]
            if i < len(parts) - 1:
                test_index = parts[i + 1]['path'].split(f'{parts_file[0]}/part-')[1].split('-')[0]
                if test_index == index:
                    continue
            file_sizes.append(x['size_bytes'])
    metadata_file = [x['path'] for x in rows_files if x['path'].endswith('metadata.json.gz')]
    if metadata_file:
        with hl.hadoop_open(metadata_file[0], 'rb') as f:
            rows_meta = json.loads(f.read())
            try:
                partition_bounds = [
                    (x['start']['locus']['contig'], x['start']['locus']['position'],
                     x['end']['locus']['contig'], x['end']['locus']['position'])
                    for x in rows_meta['jRangeBounds']]
            except KeyError:
                pass
    return partition_bounds, file_sizes


def pair_plot(
        data: pd.DataFrame,
        label_col: str = None,
        colors: Union[List[str], Dict[str, str]] = None,
        tools: str = "save,pan,box_zoom,reset,wheel_zoom,box_select,lasso_select,help",
        tooltip_cols: List[str] = None
) -> Column:
    """

    Plots each column of `data` against each other and returns a grid of plots.
    The diagonal contains a histogram of each column, or a density plot if labels are provided.
    The lower diagonal contains scatter plots of each column against each other.
    The upper diagonal is empty.
    All columns should be numerical with the exception of the `label_col` if provided.
    If a color dict containing provided mapping labels to specific colors can be specified using `color_dict`

    :param DataFrame data: Dataframe to plot
    :param str label_col: Column of the DataFrame containing the labels
    :param list of str or dict of str -> str colors: RGB hex colors. If a dict is provided, it should contain the mapping of label to colors.
    :param str tools: Tools for the resulting plots
    :param list of str tooltip_cols: Additional columns that should be displayed in tooltip
    :return: Grid of plots (column of rows)
    :rtype: Column
    """
    from bokeh.palettes import viridis

    if tooltip_cols is None:
        tooltip_cols = [] if label_col is None else [label_col]
    elif label_col not in tooltip_cols:
        tooltip_cols.append(label_col)

    if label_col is None and colors is not None:
        logger.warn('`colors_dict` ignored since no `label_col` specified')

    colors_col = '__pair_plot_color'

    colors_dict = {}
    if label_col is None:
        data[colors_col] = viridis(1) * len(data)
    else:
        if not isinstance(colors, dict):
            labels = set(data[label_col])
            color_palette = viridis(len(labels)) if colors is None else colors
            colors_dict = {l: color_palette[i] for i, l in enumerate(labels)}
        else:
            colors_dict = colors
        data[colors_col] = [colors_dict.get(l, 'gray') for l in data[label_col]]
        tools = 'hover,' + tools

    data_cols = [c for c in data.columns if c not in [colors_col, label_col] + tooltip_cols]
    data_ranges = [DataRange1d(start=rmin - (abs(rmin - rmax) * 0.05), end=rmax + (abs(rmin - rmax) * 0.05)) for rmin, rmax in zip(data[data_cols].min(axis=0), data[data_cols].max(axis=0))]
    data_source = ColumnDataSource(data={c: data[c] for c in data.columns})

    n_cols = len(data_cols)

    plot_grid = []
    for i in range(n_cols):
        row = [None] * n_cols
        for j in range(i + 1):
            p = figure(
                x_axis_label=data_cols[j] if i == n_cols - 1 else '',
                y_axis_label=data_cols[i] if j == 0 else '',
                tools=tools
            )
            p.x_range = data_ranges[j]

            if i == j:
                if label_col is None:
                    hist, edges = np.histogram(data[data_cols[i]], density=False, bins=50)
                    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
                else:
                    density_data = data[[label_col, data_cols[i]]].groupby(label_col).apply(lambda df: np.histogram(df[data_cols[i]], density=True, bins=20))
                    for label, (hist, edges) in density_data.iteritems():
                        line_source = ColumnDataSource({'edges': edges[:-1], 'hist': hist, label_col: [label] * len(hist)})
                        p.line('edges', 'hist', color=colors_dict[label], legend=label_col, source=line_source)
                        p.select_one(HoverTool).tooltips = [(label_col, f'@{label_col}')]
            else:
                p.y_range = data_ranges[i]
                if label_col is not None:
                    p.circle(data_cols[j], data_cols[i], source=data_source, color=colors_col, legend=label_col)
                else:
                    p.circle(data_cols[j], data_cols[i], source=data_source, color=colors_col)
                if tooltip_cols:
                    p.select_one(HoverTool).tooltips = [(x, f'@{x}') for x in tooltip_cols]

            row[j] = p
        plot_grid.append(row)

    return gridplot(plot_grid, toolbar_location='left')
