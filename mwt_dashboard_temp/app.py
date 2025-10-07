import streamlit as st
from config import pages, metric_palette, config
from utils.data_loader import fetch_data
from utils.preprocess import select_datasets
from app_pages import home, gene, allele, help, citations, custom_gene, custom_allele, psa, clustering
from utils.auth import check_password

# if not check_password():
#     st.stop()

st.set_page_config(page_title="MWT Dashboard", layout="wide")
st.title("MWT Data Dashboard - Rankin Lab @ UBC")

pages = [
    "Home - Getting Started",
    "Data at a Glance",
    "Gene-specific Data",
    "Allele-specific Data",
    "Custom Gene Selection",
    "Custom Allele Selection",
    "Post Stimulus Data", 
    "Gene Clustering",
    "Citations"
]

# Select a page
page = st.sidebar.radio("Select a page", pages)

# Load and filter data
data = fetch_data()
if page not in [pages[0], pages[8]]: # IMP: If page is not "Help" or "Citations" from config.py
    select_datasets(data)
data["metric_palette"] = metric_palette 
data["plotly_config"] = config

# Route pages
if page == pages[0]:
    help.render()
elif page == pages[1]:
    home.render(data)
elif page == pages[2]:
    gene.render(data)
elif page == pages[3]:
    allele.render(data)
elif page == pages[4]:
    custom_gene.render(data)
elif page == pages[5]:
    custom_allele.render(data)
elif page == pages[6]:
    psa.render(data)
elif page == pages[7]:
    clustering.render(data)
elif page == pages[8]:
    citations.render()