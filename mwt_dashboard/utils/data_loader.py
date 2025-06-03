# utils/fetch_data.py
import psycopg
import streamlit as st
import numpy as np
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
                        - id_data
    """
    with psycopg.connect(
        dbname="mwtdata", 
        user=st.secrets["psql_user"], 
        password=st.secrets["psql_passwword"], 
        host="rds-mwt-data.ctie02ksmcqc.ca-central-1.rds.amazonaws.com", 
        port=5432
        ) as connection:
        
        # Read data from SQLite database
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