
import hail as hl
import numpy as np
from ipywidgets import interact
import math
import pandas as pd
from collections import OrderedDict

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

variant_annotation_colors = {
    'synonymous_variant': COLOR_SYN,
    'missense_variant': COLOR_MIS,
    'stop_gained': COLOR_LOF,
    'Synonymous': COLOR_SYN,
    'Missense': COLOR_MIS,
    'LoF': COLOR_LOF,
}

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


def plot_hail_hist(hist_data: hl.Struct, title: str = 'Plot', log: bool = False) -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))
    """
    low = int(log)
    distance = abs(hist_data.bin_edges[0] - hist_data.bin_edges[1])
    num_data_points = sum(hist_data.bin_freq)
    title = f'{title} ({num_data_points:,} data points)'
    p = figure(title=title, y_axis_type="log") if log else figure(title=title)
    hist_data.bin_freq = [x + low for x in hist_data.bin_freq]
    p.quad(top=hist_data.bin_freq, bottom=low,
           left=hist_data.bin_edges[:-1], right=hist_data.bin_edges[1:], fill_color="#036564", line_color="#033649")
    p.quad(top=[hist_data.n_smaller + low, hist_data.n_larger + low], bottom=[low, low],
           left=[hist_data.bin_edges[0] - distance, hist_data.bin_edges[-1]],
           right=[hist_data.bin_edges[0], hist_data.bin_edges[-1] + distance])
    return p


def plot_hail_hist_cumulative(hist_data: hl.Struct, title: str = 'Plot', normalize: bool = True,
                              line_color: str = "#036564", line_width: int = 3, log: bool = False) -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))
    """
    cumulative_data = np.cumsum(hist_data.bin_freq) + hist_data.n_smaller
    np.append(cumulative_data, [cumulative_data[-1] + hist_data.n_larger])
    num_data_points = max(cumulative_data)

    if normalize: cumulative_data = cumulative_data / num_data_points
    title = f'{title} ({num_data_points:,} data points)'
    p = figure(title=title, y_axis_type="log") if log else figure(title=title)
    p.line(hist_data.bin_edges[:-1], cumulative_data, line_color=line_color, line_width=line_width)
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


def pair_plot(
        data: pd.DataFrame,
        label_col: str = None,
        colors_dict: Dict[str, str] = None,
        tools: str = "save,pan,box_zoom,reset,wheel_zoom,box_select,lasso_select,help"
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
    :param dict of str -> str colors_dict: Mapping of label to colors
    :param str tools: Tools for the resulting plots
    :return: Grid of plots (column of rows)
    :rtype: Column
    """

    if label_col is None and colors_dict is not None:
        logger.warn('`colors_dict` ignored since no `label_col` specified')

    colors_col = '__pair_plot_color'

    if label_col is None:
        data[colors_col] = ['#1f77b4'] * len(data)
    else:
        if colors_dict is None:
            from bokeh.palettes import Spectral6
            colors_dict = {l: Spectral6[i] for i,l in enumerate(set(data[label_col]))}
        data[colors_col] = [colors_dict.get(l, 'gray') for l in data[label_col]]

    data_cols = [c for c in data.columns if c not in [colors_col, label_col]]
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
                    density_data = data[[colors_col, data_cols[i]]].groupby(colors_col).apply(lambda df: np.histogram(df[data_cols[i]], density=True, bins=20))
                    for color, (hist, edges) in density_data.iteritems():
                        p.line(edges[:-1], hist, color=color)
            else:
                p.y_range = data_ranges[i]
                if label_col is not None:
                    p.circle(data.columns[j], data.columns[i], source=data_source, color=colors_col, legend=label_col)
                else:
                    p.circle(data.columns[j], data.columns[i], source=data_source, color=colors_col)

            row[j] = p
        plot_grid.append(row)

    return gridplot(plot_grid)
