# utils/fetch_data.py
import psycopg
import streamlit as st
import numpy as np
from utils.helpers import read, aggregate_unique_values, aggregate_unique_values_MSD



def fetch_data():
    """
    Connects to PostgreSQL database and loads all necessary tables for the dashboard.
    
    Returns:
        pd.DataFrame: Contains:
                        - tap_output
                        - tap_tstat_allele
                        - tap_tstat_data
                        - gene_profile_data
                        - allele_profile_data
                        - gene_MSD
                        - allele_MSD
                        - id_data
    """
    with psycopg.connect(
        dbname="mwtdata", 
        user=st.secrets["psql_user"], 
        password=st.secrets["psql_passwword"], 
        host="rds-mwt-data.ctie02ksmcqc.ca-central-1.rds.amazonaws.com", 
        port=5432
        ) as connection:
        
        data = {
            "tap_output": read('tap_response_data', connection),
            "tap_tstat_allele": aggregate_unique_values(read('tstat_gene_data', connection),["Gene"]).explode('Screen').reset_index(drop=True),
            "tap_tstat_data": aggregate_unique_values(read('tstat_allele_data', connection),["dataset"]).explode('Screen').reset_index(drop=True),
            "gene_profile_data": aggregate_unique_values(read('gene_profile_data', connection),['Gene','Metric']).explode('Screen').reset_index(drop=True),
            "allele_profile_data": aggregate_unique_values(read('allele_profile_data', connection),['dataset','Metric']).explode('Screen').reset_index(drop=True),
            "gene_MSD": aggregate_unique_values_MSD(read('gene_MSD', connection),["Gene"]).explode('Screen').reset_index(drop=True),
            "allele_MSD": aggregate_unique_values_MSD(read('allele_MSD', connection),["dataset"]).explode('Screen').reset_index(drop=True),
            "id_data": read('Gene_Allele_WormBaseID', connection) ##table in database with wormbase id's for all genes and alleles
        }
        
        data["tap_output"]["Strain"] = data["tap_output"]["Gene"] + " (" + data["tap_output"]["Allele"] + ")"
        
    return data



def select_datasets(data):
    """
    Allows user to select datasets to filter across multiple global dataframes.
    Updates global dataframes by filtering on user selection and cleaning data.

    Globals modified:
        gene_MSD, phenotype_list, tap_output, tap_tstat_allele,
        tap_tstat_data, gene_profile_data, allele_profile_data, allele_MSD
    """

    datasets = st.multiselect(
        label="Select Datasets",
        options=data["gene_MSD"].Screen.unique(),
        default="PD_Screen",
        placeholder="make a selection",
        help="select and de-select datasets you want to analyze",
        key="datasetselection"
    )

    data["datasets"] = datasets

    # Construct phenotype list
    phenotype_list = []
    for col in data["gene_MSD"].columns[1:]:
        col_split = col.split("-", 1)[0]
        phenotype_list.append(col_split)
    phenotype_list.remove('Screen')

    dropna_features = list(np.unique(phenotype_list))
    dropna_features.remove('Spontaneous Recovery of Response Duration')
    dropna_features.remove('Spontaneous Recovery of Response Probability')
    dropna_features.remove('Spontaneous Recovery of Response Speed')
    dropna_features.remove('Memory Retention of Response Duration')
    dropna_features.remove('Memory Retention of Response Probability')
    dropna_features.remove('Memory Retention of Response Speed')

    data["phenotype_list"] = phenotype_list

    # Apply filtering and cleaning
    data["tap_output"] = data["tap_output"][data["tap_output"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

    data["tap_tstat_allele"] = (
        data["tap_tstat_allele"][data["tap_tstat_allele"]["Screen"].isin(datasets)]
        .dropna(subset=dropna_features)
        .drop(columns=["Screen"])
        .replace(["N2_N2", "N2_XJ1"], "N2")
    )

    data["tap_tstat_data"] = (
        data["tap_tstat_data"][data["tap_tstat_data"]["Screen"].isin(datasets)]
        .dropna(subset=dropna_features)
        .drop(columns=["Screen"])
        .replace(["N2_N2", "N2_XJ1"], "N2")
    )

    data["gene_profile_data"] = data["gene_profile_data"][data["gene_profile_data"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["allele_profile_data"] = data["allele_profile_data"][data["allele_profile_data"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["gene_MSD"] = data["gene_MSD"][data["gene_MSD"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
    
    data["allele_MSD"] = data["allele_MSD"][data["allele_MSD"]["Screen"].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

    return data

