import streamlit as st

def render(data=None):
    st.markdown("""
    ## Home - Getting Started

    ### About the Data Dashboard:
    Since the inception of the [Multi-Worm Tracker (MWT)](https://doi.org/10.1038/nmeth.1625), the Rankin Lab has collected vast amounts of data from a large number of strains across several experimental screens. 
    This dashboard, written by Joseph Liang (PhD Candidate), consolidates MWT data collected over the years following the standard 10-second ISI (30 stimuli at 10-second inter-stimulus intervals) habituation protocols.

    Although the MWT is powerful in capturing rich multi-phenotype data from large populations of living animals simultaneously, the raw output is difficult to access for users unfamiliar with its inner workings.
    This tool addresses that barrier by processing and visualizing the appropriate datasets through an easy-to-use dashboard.

    ### The Data
    Data in this dashboard comes from a variety of experiments conducted by students of the Rankin Lab. These studies examine *C. elegans* strains with (mainly loss-of-function) mutations in genes of interest.

    All experiments use the 10s-ISI tap-habituation protocol. Heatmaps and phenomic profiles are calculated by comparing mutant strains to N2 wild-type controls run on the same days.

    ### Getting Started
    Navigate between pages using the sidebar. Gene- and allele-specific pages display summarized data for each gene or allele.
    To make cross-gene comparisons, use the “Custom Gene/Allele Selection” pages.

    You can download all plots and datasets directly from the interface. Don’t like our visualizations? Export the data and create your own using your favorite tools.

    ### Contact
    For bugs, issues, or feature requests, please contact Joseph Liang at **joseph.liang@psych.ubc.ca**.
    """)
