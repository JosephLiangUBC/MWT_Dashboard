import streamlit as st
from config import pages, metric_palette, config
from utils.data_loader import fetch_data
from utils.preprocess import select_datasets
from app_pages import home, gene, allele, help, citations, custom_gene, custom_allele, psa
from utils.auth import check_password

# if not check_password():
#     st.stop()

st.set_page_config(page_title="MWT Dashboard", layout="wide")
st.title("MWT Data Dashboard - Rankin Lab @ UBC")

# Select a page
page = st.sidebar.radio("Select a page", pages)

# Load and filter data
data = fetch_data()
if page not in [pages[6], pages[7]]: # IMP: If page is not "Help" or "Citations" from config.py
    select_datasets(data)
data["metric_palette"] = metric_palette 
data["plotly_config"] = config

# Route pages
if page == pages[0]:
    home.render(data)
elif page == pages[1]:
    gene.render(data)
elif page == pages[2]:
    allele.render(data)
elif page == pages[3]:
    custom_gene.render(data)
elif page == pages[4]:
    custom_allele.render(data)
elif page == pages[5]:
    psa.render(data)
elif page == pages[6]:
    help.render()
elif page == pages[7]:
    citations.render()