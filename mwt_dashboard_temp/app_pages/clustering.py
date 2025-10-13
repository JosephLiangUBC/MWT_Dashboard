# pages/clustering.py
import streamlit as st
import io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# clustering specific imports
from sklearn.compose import make_column_transformer
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import plotly.express as px


    
def render(data):

    st.header("Clustering Genes")

    n_clusters = st.slider(
        "Select number of clusters",
        min_value=2,
        max_value=8,
        value=6,  # default
        step=1
        )
    
    # Load data
    df = data['tap_tstat_data']

    gene_names = df['Gene'].values

    df = df[~df.isin([np.inf, -np.inf]).any(axis=1)]

    # Preprocess numerical columns 
    numeric_cols = df.select_dtypes(include='number').columns
    si = SimpleImputer(strategy="constant", fill_value=-9999)
    scaler = StandardScaler()
    numerical_transformer = make_pipeline(si, scaler)

    # Preprocess categorical columns
    categorical_cols = ['Gene']
    ohe = OneHotEncoder(handle_unknown="ignore", sparse_output=False, dtype = int)
    categorical_transformer = make_pipeline(ohe)

    # Apply column transformations
    preprocessor = make_column_transformer(
        (numerical_transformer, numeric_cols),
        (categorical_transformer, categorical_cols),
        remainder="passthrough" 
    )

    X = preprocessor.fit_transform(df)

    # PCA + KMeans
    pca = PCA(n_components=99)
    pca_df = pd.DataFrame(pca.fit_transform(X))
    kmeans = KMeans(n_clusters=n_clusters, n_init='auto', random_state=100)
    labels = kmeans.fit_predict(pca_df)

    # tSNE to visualise multidimensional data in 2D
    tsne = TSNE(n_components=2, random_state=100, perplexity=30)
    tsne_df = pd.DataFrame(tsne.fit_transform(pca_df), columns=['tSNE1', 'tSNE2'])
    tsne_df['Cluster'] = labels.astype(str)
    tsne_df['Gene'] = gene_names

    # internal tests 
    # st.write(df.head(1))
    # st.write(pca_df.head(1))
    # st.write(tsne_df.head(1))
    # st.write(tsne_df['Gene'].isnull().sum())

    # Plot 
    fig = px.scatter(
        tsne_df,
        x="tSNE1",
        y="tSNE2",
        color="Cluster",
        hover_data={"Gene": True, "Cluster": True, "tSNE1":False, "tSNE2":False},
        title="t-SNE Visualization of Gene Clusters (KMeans labels)",
        color_discrete_sequence=px.colors.qualitative.Plotly
    )

    fig.update_traces(
        marker=dict(size=6, opacity=0.8, line=dict(width=0))
        )
    
    fig.update_layout(
        plot_bgcolor='white',      
        xaxis=dict(showticklabels=False, showgrid=False, gridcolor='#eeeeee', gridwidth=0.5),
        yaxis=dict(showticklabels=False, showgrid=True, gridcolor='#eeeeee', gridwidth=0.5),
        height=600
    )

    st.plotly_chart(fig, use_container_width=True, config=data["plotly_config"])
