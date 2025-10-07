# pages/clustering.py
import streamlit as st
import io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from utils.helpers import convert_df, read
from config import config

# clustering specific imports
from sqlalchemy import create_engine
from sklearn.compose import make_column_transformer
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def render(data):
    st.header("Clustering data by genes")

    df = data['tap_tstat_data']

    df = df[~df.isin([np.inf, -np.inf]).any(axis=1)]

    numeric_cols = df.select_dtypes(include='number').columns
    si = SimpleImputer(strategy="constant", fill_value=-9999)
    scaler = StandardScaler()
    numerical_transformer = make_pipeline(si, scaler)

    # Categorical Columns
    categorical_cols = ['Gene']
    ohe = OneHotEncoder(handle_unknown="ignore", sparse_output=False, dtype = int)
    categorical_transformer = make_pipeline(ohe)

    # Apply transformations
    preprocessor = make_column_transformer(
        (numerical_transformer, numeric_cols),
        (categorical_transformer, categorical_cols),
        remainder="passthrough" 
    )

    X = preprocessor.fit_transform(df)

    # pca = PCA()
    # pca.fit(X)

    # PCA + KMeans
    pca = PCA(n_components=99)
    pca_df = pca.fit_transform(X)
    kmeans = KMeans(n_clusters=6, n_init='auto', random_state=100)
    labels = kmeans.fit_predict(pca_df)
    
    df['Cluster'] = labels


    # Plot
    pca_df_vis = PCA(n_components=2, random_state=100) # n_components=2 for 2D plot
    labels_vis = pca_df_vis.fit_transform(pca_df)                     
    centers_vis = pca_df_vis.transform(kmeans.cluster_centers_)
    method = "PCA(2)"

    plt.figure(figsize=(8,6))
    plt.scatter(labels_vis[:, 0], labels_vis[:, 1], c=labels, cmap="tab10", alpha=0.7, s=15, linewidths=0)
    plt.scatter(centers_vis[:, 0], centers_vis[:, 1], marker="X", s=10, c="black", label="Centroid")
    plt.title(f"Clusters (assigned in 100D)")
    # plt.xlabel("PCA 1"); plt.ylabel("PCA 2")
    plt.legend()
    plt.tight_layout()
    plt.show()






    # # ------ fix
    # st.write("Columns in df:", df.columns.tolist())
    # st.write("Numeric cols:", numeric_cols)
    # st.write("Categorical cols:", categorical_cols)
    # st.write("X shape:", X.shape)
    # st.write("X type:", type(X))
    # st.write("First row:", X[0])
