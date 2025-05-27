# utils/helpers.py

import pandas as pd
import sqlite3
import streamlit as st

def convert_df(df):
    """
    Convert a DataFrame to CSV format for download buttons.

    Parameters:
    df (pd.DataFrame): DataFrame to convert

    Returns:
    bytes: CSV-encoded data
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


def aggregate_unique_values(df, groupby_cols):
    """
    Groups the dataframe by the specified columns and aggregates 'Screen' values as a list of unique items.
    
    Args:
        df (pd.DataFrame): Input dataframe.
        groupby_cols (list): Columns to group by.
    
    Returns:
        pd.DataFrame: Grouped dataframe with unique 'Screen' values aggregated.
    """
    return df.groupby(groupby_cols).agg({"Screen": lambda x: list(set(x))}).reset_index()

def aggregate_unique_values_MSD(df, groupby_cols):
    """
    Groups the MSD dataframe by the specified columns and aggregates all columns except groupby_cols 
    and 'Screen' as mean, and aggregates 'Screen' as a list of unique items.
    
    Args:
        df (pd.DataFrame): Input MSD dataframe.
        groupby_cols (list): Columns to group by.
    
    Returns:
        pd.DataFrame: Aggregated MSD dataframe.
    """
    value_cols = [col for col in df.columns if col not in groupby_cols + ['Screen']]
    return df.groupby(groupby_cols).agg(
        {**{col: 'mean' for col in value_cols}, 'Screen': lambda x: list(set(x))}
    ).reset_index()
