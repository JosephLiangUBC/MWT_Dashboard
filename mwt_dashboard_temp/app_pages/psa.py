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

    psa_data = data["psa_output"]

    st.subheader("Raw Data Preview")
    st.dataframe(psa_data)

    with st.expander("Show column names and data types"):
        st.write(psa_data.dtypes)

    st.download_button(
        label="Download CSV",
        data=convert_df(psa_data),
        file_name="psa_data.csv",
        mime="text/csv"
    )