# utils/fetch_data.py
import psycopg
import streamlit as st
import numpy as np
import pandas as pd
from utils.helpers import read, aggregate_unique_values, aggregate_unique_values_MSD

@st.cache_data
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
                        - psa_output
                        - id_data
    """
    with psycopg.connect(
        dbname="mwtdata", 
        user=st.secrets["psql_user"], 
        password=st.secrets["psql_passwword"], 
        host="rds-mwt-data.ctie02ksmcqc.ca-central-1.rds.amazonaws.com", 
        port=5432
        ) as connection:
        
        # ------- Read data from PostgreSQL database -----------

        # (1) Tap Response
        tap_output =  read('tap_response_data', connection)
        tap_output["Strain"] = tap_output["Gene"] + " (" + tap_output["Allele"] + ")"
        

        # (2) Tstat: Baseline + Tap + PSA tstat data by Allele
        tap_tstat_allele = aggregate_unique_values(read('tstat_gene_data', connection),["Gene"]).explode('Screen').reset_index(drop=True)
        # normalise the data
        numeric_col = tap_tstat_allele.select_dtypes(include=np.number).columns
        tap_tstat_allele[numeric_col] = (tap_tstat_allele[numeric_col] - tap_tstat_allele[numeric_col].mean()) / tap_tstat_allele[numeric_col].std()
        # if "Neuron_Genes_Screen" in tap_tstat_allele["Screen"].values:
        #     control_row = tap_tstat_allele.loc["N2_XJ1"]
        # else:
        #     control_row = tap_tstat_allele.loc["N2"]
        #     tap_tstat_allele = tap_tstat_allele.subtract(control_row, axis=1)
        

        # (3) Tstat: Baseline + Tap + PSA tstat data by Gene
        tap_tstat_data = aggregate_unique_values(read('tstat_allele_data', connection),["dataset"]).explode('Screen').reset_index(drop=True)
        # normalise the data
        numeric_cols = tap_tstat_data.select_dtypes(include=np.number).columns
        tap_tstat_data[numeric_cols] = (tap_tstat_data[numeric_cols] - tap_tstat_data[numeric_cols].mean()) / tap_tstat_data[numeric_cols].std()
        # tap_tstat_data = tap_tstat_data.subtract(tap_tstat_data.loc["N2"], axis=1)


        # (4) MSD: Baseline + Tap + PSA by Gene
        gene_MSD = aggregate_unique_values_MSD(read('gene_MSD', connection),["Gene"]).explode('Screen').reset_index(drop=True)
        
        # (5) MSD: Baseline + Tap + PSA by Allele
        allele_MSD = aggregate_unique_values_MSD(read('allele_MSD', connection),["dataset"]).explode('Screen').reset_index(drop=True)
        
        
        # (6) Allele Profile (tstat melted) 
        gene_profile_data=tap_tstat_data.reset_index()
        gene_profile_data=pd.melt(gene_profile_data, id_vars=["dataset", "Screen"],
                                    var_name='Metric',
                                    value_name='T_score')
        
        
        # (7) Gene Profile (tstat melted) 
        allele_profile_data=tap_tstat_allele.reset_index()
        allele_profile_data=pd.melt(allele_profile_data, id_vars=["Gene", "Screen"],
                                    var_name='Metric',
                                    value_name='T_score')
        
        
        # (8) PSA summarised data
        psa_output =  read('psa_summarised_data', connection)
        
        # (9) ID data
        id_data = read('Gene_Allele_WormBaseID', connection) ##table in database with wormbase id's for all genes and alleles

        # Melted/Profile data to be read after normalisation
        # "gene_profile_data": aggregate_unique_values(read('gene_profile_data', connection),['Gene','Metric']).explode('Screen').reset_index(drop=True),
        # "allele_profile_data": aggregate_unique_values(read('allele_profile_data', connection),['dataset','Metric']).explode('Screen').reset_index(drop=True),
        

        # ------------- Package the datasets --------------

        data = {
            "tap_output": tap_output,
            "psa_output": psa_output,
            "tap_tstat_allele": tap_tstat_allele,
            "tap_tstat_data": tap_tstat_data,
            "gene_profile_data": gene_profile_data,
            "allele_profile_data": allele_profile_data,
            "gene_MSD": gene_MSD,
            "allele_MSD": allele_MSD,
            "id_data": id_data
        }


    return data