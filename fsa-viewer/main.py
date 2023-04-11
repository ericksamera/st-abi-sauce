#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for sanger-sequence-trim.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "stable enough"
# --------------------------------------------------
import streamlit as st
import streamlit_ext as ste
import pandas as pd
from pathlib import Path
# --------------------------------------------------
import plotly.graph_objects as go
import math
from Bio import SeqIO
from datetime import datetime
from io import StringIO
import zipfile
import io
import copy
import time
# --------------------------------------------------
class App:
    def __init__(self):
        """
        """
        self._loaded_abi_files = {}
    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        title = "fsa-viewer"
        st.set_page_config(
            page_title=f"abi-sauce | {title}",
            page_icon=':apple:',
            layout='wide',
            initial_sidebar_state='collapsed')
        st.title(title)
        st.markdown('This script is intended for processing a `.ab1` files into Mott algorithm-trimmed FASTAs.')
        st.markdown('Check out the better-maintained command-line interface on [GitHub](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger-sequence-trim.py)!')
        self._init_file_uploader()
        self._init_sidebar()
    def _init_file_uploader(self) -> None:
        """
        Function instantiates file uploading and processing for the files.
        """

        form_instance = st.empty()
        alignment_result = None

        with form_instance.form("sanger upload form", clear_on_submit=True):
            uploaded_file: list = st.file_uploader(
                'Upload `.fasta` files.',
                type=['.csv'],
                accept_multiple_files=False)

            submitted = st.form_submit_button("Align and generate consensus")

        if submitted and uploaded_file:
            st.session_state._input_csv_dict = {}
            for line in uploaded_file.getvalue().decode("utf-8").split('\n')[1:]:
                line_info = [thing.replace('"', '') for thing in line.strip().split(',')]
                if not len(line_info) > 7: continue
                color = line_info[0].lower()

                if color not in st.session_state._input_csv_dict:
                    st.session_state._input_csv_dict[color] = {
                        'x': [],
                        'y': [],
                    }

                st.session_state._input_csv_dict[color]['x'] += [float(line_info[8]), float(line_info[2]), float(line_info[10])] + [None]
                st.session_state._input_csv_dict[color]['y'] += [0, float(line_info[3]), 0] + [None]
            with st.spinner('Trimming...'):
                time.sleep(0.5)
            form_instance.empty()
            self._plot()

            if submitted and not uploaded_file:
                st.error('Attach some files first!')

    def _plot(self):
        """
        """

        max_height = max([yval for yval in st.session_state._input_csv_dict['orange']['y'] if yval]) + 200
        fig = go.Figure()
        for color in st.session_state._input_csv_dict:
            name_val = 'LADDER' if color == 'orange' else color
            color_val = color if color != 'yellow' else 'black'
            fig.add_trace(
                go.Scatter(
                    x=st.session_state._input_csv_dict[color]['x'],
                    y=st.session_state._input_csv_dict[color]['y'],
                    mode='lines',
                    #hoverinfo=hoverval,
                    name=name_val.upper(),
                    marker=dict(
                        color=color_val,
                    )
                ))
        fig.update_layout(
            height=800,
            xaxis=dict(constrain='domain'),
            yaxis=dict(constrain='range', range=[-1, max_height]),
        )
        st.plotly_chart(fig, use_container_width=True)

    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
        
        with st.sidebar:
            st.header('About')
            st.markdown(
                'In Sanger sequencing, the beginning and end of the ' 
                'electropherogram generally end up a litle messy due to the '
                'inherent randomness of the chain-termination method.')
        
            st.markdown(
                'The low-quality basecalls at the beginning and end of the '
                'electropherogram are likely not real -- and therefore not '
                'useful for BLAST or variant identification.')

            st.markdown(
                'These are usually trimmed off manually, but you can use this '
                'tool to do it too!')
            st.caption('')
            st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()