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

    psa_df = data["psa_output"]

    # Filter boxes for metric and summary
    metric_option = st.selectbox("Select metric", ['Instantaneous Speed', 'Bias', 'Angular Speed', 'Kink', 'Crab', 'Aspect Ratio', 'Curve'])
    summary_option = st.selectbox("Select summary", ['Initial', 'Recovery', 'Peak', 'Peak Tap Number', 'Initial_to_peak', 'Peak_to_recovery', 'Average'])

    # # Filter box for gene
    # gene_option = st.selectbox("Select a gene", psa_df["Gene"].unique())
    # psa_df = psa_df[psa_df["Gene"] == gene_option]

    # sort by gene
    gene_order = (
    psa_df.groupby("Gene")[f"{summary_option} PSA {metric_option}"]
    .mean()
    .sort_values()
    .index
    )

    # Plot
    st.subheader(f"Full Post-Stimulus Arousal Response — {metric_option}")
    
    fig, ax = plt.subplots(figsize=(15, 5))
    sns.set_context("poster")
    sns.barplot(
        x="Gene",
        y=f"{summary_option} PSA {metric_option}",
        data=psa_df,
        hue = "Gene",
        order = gene_order,
        estimator='mean',
        errorbar=("ci", 95),
        ax=ax,
        orient="v",
        legend=False
    )
    # sns.stripplot(
    #     x="Tap_num",
    #     y=f"{metric_option}",
    #     data=df_filtered,
    #     size=6,
    #     color="k",
    #     ax=ax
    # )
    ax.set_title(f"Post Stimulus Arousal {summary_option} {metric_option} - all genes", fontsize=14)
    ax.set_ylabel(f"{summary_option} {metric_option}", fontsize=12)
    ax.set_xlabel("Gene", fontsize=12)
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=12)
    # ax.set_ylim(top=0.33)
    st.pyplot(fig)







    # Summary charts?
    # st.subheader(f"{metric_option} by Summary Type")

    # summaries = ['Initial', 'Recovery', 'Peak', 'Initial_to_peak', 'Peak_to_recovery', 'Average']
    # col1, col2 = st.columns(2)

    # for i, summary in enumerate(summaries):
    #     col = col1 if i % 2 == 0 else col2

    #     # Build the correct column name
    #     column_name = f"PSA {summary} {metric_option}"
    #     if column_name not in df_filtered.columns:
    #         col.warning(f"{column_name} not found in data.")
    #         continue

    #     fig, ax = plt.subplots(figsize=(7, 5))
    #     sns.set_context("talk")
    #     sns.barplot(
    #         x="Screen", y=column_name,
    #         data=df_filtered,
    #         hue="Screen",
    #         ax=ax,
    #         estimator='mean',
    #         errorbar=("ci", 95),
    #         palette="muted",
    #         legend=False
    #     )
    #     sns.stripplot(
    #         x="Screen", y=column_name,
    #         data=df_filtered,
    #         size=5,
    #         color="k",
    #         ax=ax
    #     )
    #     ax.set_title(f"{summary} {metric_option}")
    #     ax.set_xlabel("")
    #     ax.set_ylabel(metric_option)
    #     ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    #     col.pyplot(fig)