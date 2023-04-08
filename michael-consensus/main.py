#!/usr/bin/env python3
__description__ =\
"""
Purpose: Streamlit wrapper for Michael's consensus implementation.
"""
__author__ = "Erick Samera"
__version__ = "2.5.0"
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
            uploaded_files: list = st.file_uploader(
                'Upload `.fasta` files.',
                type=['.fasta'],
                accept_multiple_files=True)

            submitted = st.form_submit_button("Align and generate consensus")
            if submitted:
                fasta_entries = []
                for file in uploaded_files:
                    for fasta_entry in SeqIO.parse(StringIO(file.getvalue().decode("utf-8")), 'fasta'):
                        fasta_entries.append(fasta_entry)
                SeqIO.write(fasta_entries, st.session_state.temp_dir.joinpath('temp-concat.fasta'), 'fasta')
                with st.spinner('Generating alignment ...'):
                    alignment_result = subprocess.call(
                        f"mafft --auto --quiet {st.session_state.temp_dir.joinpath('temp-concat.fasta')} > {st.session_state.temp_dir.joinpath('temp-aligned.fasta')}", 
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
    def _init_plot_data(self):
        """
        """

        def _generate_consensus(alignment_object) -> str:
            """

            """

            # TODO: add some support for degenerate nucleotides in the consensus?
            consensus_str: str = ''
            left_right_boundaries: list = [(list(re.finditer('\w', str(SeqObj.seq).upper()))[0].start(), list(re.finditer('\w', str(SeqObj.seq).upper()))[-1].end()) for SeqObj in alignment_object]
            
            for pos in range(alignment_object.get_alignment_length()):
                nuc_at_pos = [seq.seq[pos].upper() for i, seq in enumerate(alignment_object) if left_right_boundaries[i][0] <= pos <= left_right_boundaries[i][1]]
                nuc_counts = Counter(nuc_at_pos)
                if len(nuc_counts) > 1:
                    counts_at_pos = nuc_counts.most_common(5)
                    primary_nuc, secondary_nuc = counts_at_pos[0][0], counts_at_pos[1][0]

                    nuc_to_add = primary_nuc[0]

                    if primary_nuc == '-':
                        if (counts_at_pos[0][1]/len(nuc_at_pos)) > st.session_state.gap_threshold/100: nuc_to_add = '-'
                        else: nuc_to_add = secondary_nuc
                    consensus_str += nuc_to_add
                else:
                    consensus_str += nuc_counts.most_common(1)[0][0]
            return consensus_str

        nucleotide_transform = {nucleotide: i for i, nucleotide in enumerate('-ACGTUWSMKRYBDHVN')}

        if 'alignment_values' not in st.session_state:
            st.session_state.fasta_iterator = [seq_object for seq_object in SeqIO.parse(st.session_state.temp_dir.joinpath('temp-aligned.fasta'),'fasta')]

            st.session_state.alignment_values: list = [[nucleotide_transform.get(char.upper()) for char in seq_object.seq] for seq_object in st.session_state.fasta_iterator]
            st.session_state.alignment_text: list = [[char for char in seq_object.seq.upper()] for seq_object in st.session_state.fasta_iterator]
            st.session_state.seq_ids: list = [seq.id for seq in st.session_state.fasta_iterator]
            st.session_state.num_seqs: int = len(st.session_state.seq_ids)

            # VALUES FOR RELATIVE PLOTTING
            st.session_state.consensus_height = 50
            st.session_state.alignment_height = 100 + (15*len(st.session_state.alignment_text))

            degenerate_nucleotides_len = max(set(sum(st.session_state.alignment_values, []))) - 5 + 1
            st.session_state.nuc_colors = ['#ffffff']+['#64f73f', '#ffb340', '#eb413c', '#3c88ee'] + ['#a800ff']*degenerate_nucleotides_len

            st.session_state.consensus_object = AlignIO.read(st.session_state.temp_dir.joinpath('temp-aligned.fasta'),'fasta')

        st.session_state.consensus_str: str = _generate_consensus(st.session_state.consensus_object)
        st.session_state.consensus_vals = [[nucleotide_transform.get(char.upper()) for char in st.session_state.consensus_str]]
        st.session_state.consensus_chars = [[char for char in st.session_state.consensus_str]]
        st.session_state._alignment_str_to_file = '\n'.join([seq.format('fasta') for seq in st.session_state.fasta_iterator])
        self._plot()

    def _plot(self):
        st.success(f'Successfully generated a multiple sequence alignment from {st.session_state.num_seqs} sequences!', icon=None)

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                            vertical_spacing = 0,
                            row_heights=[st.session_state.consensus_height, st.session_state.alignment_height]
                            )
        fig.add_trace(go.Heatmap(
                z=st.session_state.consensus_vals,
                zmin=0,
                zmax=17,
                x=list(range(len(st.session_state.consensus_chars[0]))),
                text=st.session_state.alignment_text,
                hoverinfo='skip',
                coloraxis = "coloraxis"),
                row=1, col=1)
        fig.add_trace(go.Heatmap(
                z=st.session_state.alignment_values,
                zmin=0,
                zmax=17,
                x=list(range(len(st.session_state.alignment_values[0]))),
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

        middle = int(len((st.session_state.alignment_values[0]))/2)

        fig.update_layout(
            height=st.session_state.consensus_height + (st.session_state.alignment_height*2),
            xaxis2=dict(range=[middle, middle+200], rangeslider=dict(visible=True, thickness=0.1)), 
            xaxis=dict(), 
            yaxis=dict(zeroline=False, showgrid=True, range=[-2, 2], scaleratio=1, fixedrange=True, tickvals=[0], ticktext=['Consensus (gapped)']),
            yaxis2=dict(zeroline=False, showgrid=True, scaleratio=1, fixedrange=True, tickvals=[0+k for k in range(st.session_state.num_seqs)], ticktext=st.session_state.seq_ids),
            coloraxis={'colorscale': st.session_state.nuc_colors},
            showlegend=True,
            legend=dict(itemclick=False))
        fig.update_coloraxes(showscale=False)

        with st.expander('Show consensus sequence (gapped) in FASTA format', expanded=False):
            st.code(
                SeqRecord(
                    Seq(''.join(st.session_state.consensus_chars[0])),
                    id='Consensus (gapped)',
                    description='',).format('fasta'),
                line_numbers=True)
        with st.expander('Show consensus sequence in FASTA format', expanded=False):
            st.code(
                SeqRecord(
                    Seq(''.join(st.session_state.consensus_chars[0]).replace('-', '')),
                    id='Consensus',
                    description='',).format('fasta'),
                line_numbers=True)
        
        ste.download_button(f'Download multiple-sequence alignment (MSA)', st.session_state._alignment_str_to_file, file_name='alignment.fasta')
        st.plotly_chart(fig, use_container_width=True)
# --------------------------------------------------
if __name__=="__main__":
    ct = App()
    ct._init_page()