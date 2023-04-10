#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for assembling .ab1 files.
"""
__author__ = "Erick Samera"
__version__ = "0.1.0"
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
class App:
    def __init__(self):
        """
        """

    def _init_page(self) -> None:
        """
        Function instantiates the main page.
        """
        title = "Michael's consensus"
        st.set_page_config(
            page_title=f"abi-sauce | {title}",
            page_icon=':apple:',
            initial_sidebar_state='collapsed',
            layout='wide',
            menu_items={
                'About': "Let Erick know if this is broken!"
                }
            )
        if 'alignment_values' not in st.session_state:
            st.title(f":apple: abi-sauce | {title}")
            st.markdown("Generate a `multiple sequence alignment (MSA)` and produce a consensus sequence (using Michael's preferred calculation).")
            st.session_state.gap_threshold = 65
            self._init_file_management()
            self._init_file_uploader()
            self._init_sidebar()
    def _init_file_management(self) -> None:
        """
        Function instantiates the temp dir for temporary storage of sequences
        """
        st.session_state.temp_dir: Path = Path(__file__).parent.joinpath('intermediate-files')
        if st.session_state.temp_dir.exists(): shutil.rmtree(st.session_state.temp_dir)
        st.session_state.temp_dir.mkdir(parents=True, exist_ok=True)
    def _init_file_uploader(self) -> None:
        """
        Function instantiates file uploading and processing for the files.
        """

        form_instance = st.empty()
        alignment_result = None

        with form_instance.form("sanger upload form", clear_on_submit=True):
            col1, col2, col3 = st.columns(3)
            with col1:
                stringency_selection  =st.selectbox(
                    'Trimming stringency', 
                    ['0: No trimming', '1: LEAST', '2', '3', '4', '5', '6', '7', '8', '9: MOST',],
                    index=4)
                pratio_var = st.slider('Peak ratio (to call primary and secondary bases)', 0, 100, 33, 5, help='HelpText')
            with col3:
                upload_reference: list = st.file_uploader(
                    'To use reference-based assembly, upload a reference `.fasta`.',
                    type=['.fasta'],
                    accept_multiple_files=False)

            uploaded_files: list = st.file_uploader(
                'Upload `.ab1` files to assemble.',
                type=['.ab1'],
                accept_multiple_files=True)

            submitted = st.form_submit_button("Generate contig")
            if submitted:
                fasta_entries = []
                for abi_file in uploaded_files:
                    with open(st.session_state.temp_dir.joinpath(abi_file.name), mode='wb') as temp_write:
                        temp_write.write(abi_file.read())
                with st.spinner('Generating contig ...'):
                    alignment_result = subprocess.call(
                        f"tracy assemble --trim {stringency_selection.split(':')[0]} --format fastq data/*.ab1", 
                        shell=True)
        if alignment_result == 0:
            form_instance.empty()
            self._init_plot_data()
    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
        
        with st.sidebar:
            st.header('About')
            st.markdown(
                "In Jalview's behaviour for generating a consensus sequence, "+
                "it ignores gaps as contributors to the consensus sequence. " +
                "This means that only non-gap nucleotides are used in the "   +
                "calculation for 'majority'."
                )
            st.markdown(
                "Michael prefers that gaps be counted. It makes sense that "  +
                "if the majority of sequences have a gap there, it should be "+
                "considered a gap in the consensus." 
                )
            st.caption('')
            st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')

# --------------------------------------------------
if __name__=="__main__":
    streamlit_app = App()
    streamlit_app._init_page()