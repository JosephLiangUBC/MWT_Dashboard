# utils/helpers.py

import pandas as pd
import sqlite3
import streamlit as st
import numpy as np

@st.cache_data
def convert_df(df):
    """
    Converts a DataFrame to a UTF-8 encoded CSV for downloading.
    """
    return df.to_csv(index=False).encode('utf-8')


def read(table, connection):
    """
    Fetches all rows from a specified table in a PostgreSQL database.
    
    Inputs:
        table (str): table name in PostgreSQL database
        connection (psycopg.Connection): Active psycopg database connection
    
    Returns:  
        pd.DataFrame: Table data as a DataFrame
    """
    with connection.cursor() as cursor:
        cursor.execute(f'SELECT * from "{table}"')

        # Fetch all rows from database
        record = cursor.fetchall()
        column_names = [desc[0] for desc in cursor.description]

    return pd.DataFrame(data=record, columns=column_names)


@st.cache_data
def aggregate_unique_values(df, by):
    """
    Aggregates by specified keys and computes mean of numerical columns. 
    Transforms tstat_gene_data, tstat_allele_data , gene_profile_data, 
    and allele_profile_data tables
    
    Inputs:
        df (pd.DataFrame): Input DataFrame containing data to aggregate.
        by (list of str): Columns to group by
    
    Returns:
        grouped (pd.DataFrame): Aggregated dataframe with mean values and 
                                grouped 'Screen' list
    """
    if len(by) > 1 and 'Metric' in df.columns:
        df['Metric'] = pd.Categorical(df['Metric'], categories=df['Metric'].unique(), ordered=True)

    # Perform groupby aggregation
    grouped = df.groupby(by).agg(
        {col: 'mean' for col in df.columns if col not in by + ['Screen']}
    ).reset_index()

    # Aggregate the Screen column into a list
    grouped['Screen'] = df.groupby(by)['Screen'].apply(lambda x: list(set(x))).reset_index(drop=True)

    return grouped


@st.cache_data
def aggregate_unique_values_MSD(df, by):
    """
    Aggregates dataframe with weighted means and calculates Confidence Intervals.
    Transforms gene_MSD and allele_MSD tables

    Inputs:
        df (pd.DataFrame): Input DataFrame containing containing metrics 
                           with '-mean', '-sem', and '-count' suffixes
        by (list of str): Columns to group by

    Returns:
        grouped pd.DataFrame: Aggregated DataFrame with weighted means, updated 95% CI columns, 
        and grouped 'Screen'
    """

    # Define the columns to aggregate
    agg_cols = [col for col in df.columns if col not in by + ['Screen']]

    # Aggregate the columns using a weighted average
    grouped = df.groupby(by).apply(lambda x: pd.Series({
        col: (x[col] * x[col.replace('-mean', '-count')]).sum() / x[col.replace('-mean', '-count')].sum()
        if '-mean' in col
        else np.sqrt((x[col] ** 2 * x[col.replace('-sem', '-count')]).sum() / x[col.replace('-sem', '-count')].sum())
        if '-sem' in col
        else x[col].sum()
        if '-count' in col
        else np.nan
        for col in agg_cols
    })).reset_index()

    # Calculate new confidence intervals
    for col in grouped.columns:
        if '-mean' in col:
            mean_col = col
            sem_col = col.replace('-mean', '-sem')
            count_col = col.replace('-mean', '-count')
            ci95_lo_col = col.replace('-mean', '-ci95_lo')
            ci95_hi_col = col.replace('-mean', '-ci95_hi')

            grouped[ci95_lo_col] = grouped[mean_col] - 1.96 * grouped[sem_col] / np.sqrt(grouped[count_col])
            grouped[ci95_hi_col] = grouped[mean_col] + 1.96 * grouped[sem_col] / np.sqrt(grouped[count_col])

    # Aggregate the Screen column into a list
    grouped['Screen'] = df.groupby(by)['Screen'].apply(lambda x: list(set(x))).reset_index(drop=True)

    return grouped
