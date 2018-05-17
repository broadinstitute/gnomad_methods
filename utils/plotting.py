
import hail as hl
import numpy as np
from ipywidgets import interact

import bokeh
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file
from bokeh.io import output_notebook, push_notebook, export_png
from bokeh.models.widgets import Tabs, Panel
from bokeh.models import (
    Label, LabelSet,
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar,
)

TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"


def plot_hail_hist(hist_data: hl.Struct, title: str) -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))
    """
    distance = abs(hist_data.bin_edges[0] - hist_data.bin_edges[1])
    p = figure(title=title)
    p.quad(top=hist_data.bin_freq, bottom=0, left=hist_data.bin_edges[:-1], right=hist_data.bin_edges[1:], fill_color="#036564", line_color="#033649")
    p.quad(top=[hist_data.n_smaller, hist_data.n_larger], bottom=[0, 0], left=[hist_data.bin_edges[0] - distance, hist_data.bin_edges[-1]], right=[hist_data.bin_edges[0], hist_data.bin_edges[-1] + distance])
    return p


def plot_hail_hist_cumulative(hist_data: hl.Struct, title: str, normalize: bool = True,
                              line_color: str = "#036564", line_width: int = 3) -> bokeh.plotting.Figure:
    """
    hist_data can (and should) come straight from ht.aggregate(hl.agg.hist(ht.data, start, end, bins))
    """
    cumulative_data = np.cumsum(hist_data.bin_freq) + hist_data.n_smaller
    np.append(cumulative_data, [cumulative_data[-1] + hist_data.n_larger])
    num_data_points = max(cumulative_data)

    if normalize: cumulative_data = cumulative_data / num_data_points
    p = figure(title=f'{title} ({num_data_points} data points)')
    p.line(hist_data.bin_edges[:-1], cumulative_data, line_color=line_color, line_width=line_width)
    return p
