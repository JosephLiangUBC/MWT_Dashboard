# pages/psa.py
import streamlit as st
import io
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from utils.helpers import convert_df, read
from config import config


def render(data):
    st.header("Post Stimulus Arousal Data")

    st.session_state.setdefault('allele_select', None)
    allele_option = st.selectbox(
        'Select an allele',
        [allele for allele in data["tap_output"]['dataset'].unique() if allele != 'N2'],
        key="alleleselect",
        index=[allele for allele in data["tap_output"]['dataset'].unique() if allele != 'N2'].index(st.session_state.allele_select) if st.session_state.allele_select else 0
    )
    st.session_state.allele_select = allele_option