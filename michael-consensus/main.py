#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for Michael's consensus implementation.
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
            initial_sidebar_state='collapsed',
            layout='wide')
        st.title(f":apple: abi-sauce | {title}")
        st.markdown("Generate a `multiple sequence alignment (MSA)` and produce a consensus sequence (using Michael's preferred calculation).")
        self._init_file_uploader()
        self._init_sidebar()

    def _init_file_uploader(self) -> None:
        """
        Function instantiates file uploading and processing for the files.
        """

        form_instance = st.empty()
        x = None
        with form_instance.form("sanger upload form", clear_on_submit=True):
            _uploaded_files: list = st.file_uploader(
                        'Upload `.fasta` files.',
                        type=['.fasta'],
                        accept_multiple_files=True)

            _submitted = st.form_submit_button("Align and generate consensus")
            if _submitted:
                _fasta_entries = []
                for file in _uploaded_files:
                    for fasta_entry in SeqIO.parse(StringIO(file.getvalue().decode("utf-8")), 'fasta'):
                        _fasta_entries.append(fasta_entry)
                SeqIO.write(_fasta_entries, Path(__file__).parent.joinpath('temp-concat.fasta'), 'fasta')
                x = subprocess.call(f"mafft --auto --quiet {Path(__file__).parent.joinpath('temp-concat.fasta')} > {Path(__file__).parent.joinpath('temp-aligned.fasta')}", shell=True)
        if x==0:
            form_instance.empty()
            self.plot()


    def _init_sidebar(self) -> None:
        """
        Instantiate the sidebar.
        """
        
        with st.sidebar:
            st.header('About')
            st.markdown(
                "In Jalview's behaviour for generating a consensus sequence, "+
                "it ignores gaps as contributors to the consensus sequence. "+
                "This means that only non-gap nucleotides are used in the " +
                "calculation for 'majority'."
                )
            st.markdown(
                "Michael prefers that gaps be counted. It makes sense that "+
                "if the majority of sequences have a gap there, it should be "+
                "considered a gap in the consensus." 
                )
            st.caption('')
            st.caption(f'[@{__author__}](https://github.com/ericksamera)\t|\tv{__version__}\t|\t{__comments__}')
    
    def plot(self):
        """
        """
        _translate = {nucleotide: i for i, nucleotide in enumerate('-ACGTUWSMKRYBDHVN')}

        aln_vals = [[_translate.get(char.upper()) for char in seq_object.seq] for seq_object in SeqIO.parse(Path(__file__).parent.joinpath('temp-aligned.fasta'),'fasta')]
        aln_text = [[char for char in seq_object.seq.upper()] for seq_object in SeqIO.parse(Path(__file__).parent.joinpath('temp-aligned.fasta'),'fasta')]



        seq_ids = [seq.id for seq in SeqIO.parse(Path(__file__).parent.joinpath('temp-aligned.fasta'),'fasta')]

        consensus_height = 50
        alignment_height = 100 + (15*len(aln_text))

        nucs = []
        for i in aln_vals:
            nucs += i
        degenerate_nucleotides_len = max(set(nucs)) - 5 + 1

        nuc_colors = ['#ffffff']+['#64f73f', '#ffb340', '#eb413c', '#3c88ee'] + ['#a800ff']*degenerate_nucleotides_len

        consensus_str = ''
        consensus_file = AlignIO.read(Path(__file__).parent.joinpath('temp-aligned.fasta'),'fasta')
        left_right_boundaries = [(list(re.finditer('\w', str(SeqObj.seq).upper()))[0].start(), list(re.finditer('\w', str(SeqObj.seq).upper()))[-1].end()) for SeqObj in consensus_file]

        for pos in range(consensus_file.get_alignment_length()):
            nuc_at_pos = [seq.seq[pos].upper() for i, seq in enumerate(consensus_file) if left_right_boundaries[i][0] <= pos <= left_right_boundaries[i][1]]
            if len(Counter(nuc_at_pos)) > 1:
                counts_at_pos = Counter(nuc_at_pos).most_common(5)
                primary_nuc, secondary_nuc = counts_at_pos[0][0], counts_at_pos[1][0]

                nuc_to_add = primary_nuc[0]

                if primary_nuc == '-':
                    if (counts_at_pos[0][1]/len(nuc_at_pos)) > 65/100: nuc_to_add = '-'
                    else: nuc_to_add = secondary_nuc
                consensus_str += nuc_to_add
            else:
                consensus_str += Counter(nuc_at_pos).most_common(1)[0][0]
        consensus_vals = [[_translate.get(char.upper()) for char in consensus_str]]
        consensus_chars = [[char for char in consensus_str]]

        num_sequences = len(aln_vals)

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                            vertical_spacing = 0,
                            row_heights=[consensus_height, alignment_height]
                            )

        fig.add_trace(go.Heatmap(
                z=consensus_vals,
                zmin=0,
                zmax=17,
                x=list(range(len(consensus_chars[0]))),
                text=aln_text,
                hoverinfo='skip',
                coloraxis = "coloraxis"),
                row=1, col=1)
        fig.add_trace(go.Heatmap(
                z=aln_vals,
                zmin=0,
                zmax=17,
                x=list(range(len(aln_vals[0]))),
                hoverinfo='skip',
                coloraxis = "coloraxis"),
                row=2, col=1)

        for nucleotide, trace_color in zip('ACGT', ['#64f73f', '#ffb340', '#eb413c', '#3c88ee']):
            fig.add_trace(go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    name=nucleotide,
                    marker=dict(size=20, color=trace_color, symbol='square')))

        middle = int(len((aln_vals[0]))/2)

        fig.update_layout(
            height=consensus_height + (alignment_height*2),
            xaxis2=dict(range=[middle, middle+200], rangeslider=dict(visible=True, thickness=0.1)), 
            xaxis=dict(), 
            yaxis=dict(zeroline=False, showgrid=True, range=[-2, 2], scaleratio=1, fixedrange=True, tickvals=[0], ticktext=['Consensus (gapped)']),
            yaxis2=dict(zeroline=False, showgrid=True, scaleratio=1, fixedrange=True, tickvals=[0+k for k in range(num_sequences)], ticktext=seq_ids),
            coloraxis={'colorscale':nuc_colors},
            showlegend=True,
            legend=dict(itemclick=False))
        fig.update_coloraxes(showscale=False)

        with st.expander('Show consensus sequence (gapped) in FASTA format', expanded=False):
            st.code(
                SeqRecord(
                    Seq(''.join(consensus_chars[0])),
                    id='Consensus (gapped)',
                    description='',).format('fasta'),
                line_numbers=True)
        with st.expander('Show consensus sequence in FASTA format', expanded=False):
            st.code(
                SeqRecord(
                    Seq(''.join(consensus_chars[0]).replace('-', '')),
                    id='Consensus',
                    description='',).format('fasta'),
                line_numbers=True)
        _alignment_str_to_file = '\n'.join([seq.format('fasta') for seq in SeqIO.parse(Path(__file__).parent.joinpath('temp-aligned.fasta'), 'fasta')])
        ste.download_button('Download multiple-sequence alignment (MSA)', _alignment_str_to_file, file_name='alignment.fasta')
        st.plotly_chart(fig, use_container_width=True)
        
# --------------------------------------------------
if __name__=="__main__":
    ct = App()
    ct._init_page()