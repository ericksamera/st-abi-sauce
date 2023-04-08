#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for assembling multiple Sanger electropherograms.
"""
__author__ = "Erick Samera"
__version__ = "0.0.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
from pathlib import Path
import shutil
# --------------------------------------------------
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
from io import StringIO
import subprocess
import re
# --------------------------------------------------
