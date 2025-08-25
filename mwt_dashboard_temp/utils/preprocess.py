# utils/preprocess.py
import streamlit as st
import numpy as np

def select_datasets(data):
    datasets = st.multiselect(
        label="Select Datasets",
        options=data["gene_MSD"].Screen.unique(),
        default=data["gene_MSD"].Screen.unique()[0],
        placeholder="make a selection",
        help="select and de-select datasets you want to analyze",
        key="datasetselection"
    )

    data["datasets"] = datasets

    phenotype_list = []
    for col in data["gene_MSD"].columns[1:]:
        col_split = col.split("-", 1)[0]
        phenotype_list.append(col_split)
    phenotype_list.remove('Screen')

    dropna_features = list(np.unique(phenotype_list))
    for f in [
        'Spontaneous Recovery of Response Duration',
        'Spontaneous Recovery of Response Probability',
        'Spontaneous Recovery of Response Speed',
        'Memory Retention of Response Duration',
        'Memory Retention of Response Probability',
        'Memory Retention of Response Speed'
    ]:
        dropna_features.remove(f)

    data["phenotype_list"] = phenotype_list

    data["tap_output"] = data["tap_output"][data["tap_output"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

    data["psa_output"] = data["psa_output"][data["psa_output"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

    data["tap_tstat_allele"] = (
        data["tap_tstat_allele"][data["tap_tstat_allele"]["Screen"].isin(datasets)]
        # .dropna(subset=[col for col in dropna_features if col in data["tap_tstat_allele"].columns])
        .drop(columns=["Screen"])
        .replace(["N2_N2", "N2_XJ1"], "N2")
    )

    data["tap_tstat_data"] = (
        data["tap_tstat_data"][data["tap_tstat_data"]["Screen"].isin(datasets)]
        # .dropna(subset=[col for col in dropna_features if col in data["tap_tstat_data"].columns])
        .drop(columns=["Screen", "Peak Tap Number of PSA Angular Speed"])
        .replace(["N2_N2", "N2_XJ1"], "N2")
    )

    data["gene_profile_data"] = data["gene_profile_data"][data["gene_profile_data"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["allele_profile_data"] = data["allele_profile_data"][data["allele_profile_data"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["gene_MSD"] = data["gene_MSD"][data["gene_MSD"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["allele_MSD"] = data["allele_MSD"][data["allele_MSD"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

    return data
