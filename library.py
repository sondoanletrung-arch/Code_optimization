from __future__ import annotations

import os

import numpy as np
import pandas as pd

from scipy import io

import scanpy as sc

import igraph

import infercnvpy as icv

import harmonypy

# Core scverse libraries

import anndata as ad

# Data retrieval
import pooch

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

from io import BytesIO
from PIL import Image
import base64

from bokeh.plotting import figure, show, output_notebook, save
from bokeh.models import HoverTool, ColumnDataSource, CategoricalColorMapper, Label, LabelSet
from bokeh.transform import factor_cmap
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis
from bokeh.io import export_png
import skimage 