
import hail as hl
import numpy as np
from ipywidgets import interact
import math

import bokeh
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file
from bokeh.io import output_notebook, push_notebook, export_png
from bokeh.models.widgets import Tabs, Panel
from bokeh.palettes import *
from bokeh.models import *
from typing import *
from bokeh.plotting.helpers import stack

TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

COLOR_CPG = '#2E9FFE'
COLOR_TI = '#458B00'
COLOR_TV = '#EA4444'
ACG = '#422efe'
CCG = '#2e73fe'
GCG = '#2ec3fe'
TCG = '#2efef7'

variant_type_colors = {
    'CpG': COLOR_CPG,
    'non-CpG transition': COLOR_TI,
    'transversion': COLOR_TV,
    'ACG': ACG, 'CCG': CCG, 'GCG': GCG, 'TCG': TCG
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
              'hk': '#C60C30',
              'sg': 'darkred',
              'tw': '#009246',
              'unk': COLOR_OTH}
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
             'bg': 'Bulgaria',
             'de': 'Germany',
             'ee': 'Estonia',
             'es': 'Spain',
             'fi': 'Finland',
             'gb': 'UK',
             'neu': 'Northern Europe',
             'seu': 'Southern Europe',
             'it': 'Italy',
             'se': 'Sweden',
             'cn': 'China',
             'kr': 'Korea',
             'hk': 'Hong Kong',
             'sg': 'Singapore',
             'tw': 'Taiwan',
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


def plot_hail_hist_both(hist_data: hl.Struct, title: str, normalize: bool = True,
                        line_color: str = "#036564", line_width: int = 3):
    p1 = plot_hail_hist(hist_data, title)
    p2 = plot_hail_hist_cumulative(hist_data, f'{title} (Cumulative)', normalize, line_color, line_width)
    return gridplot([[p1, p2]])


def linear_and_log_tabs(plot_func: Callable, df: Any) -> Tabs:
    panels = []
    for axis_type in ["linear", "log"]:
        fig = plot_func(df, axis_type)
        panel = Panel(child=fig, title=axis_type)
        panels.append(panel)

    return Tabs(tabs=panels)
