import streamlit as st
from config import pages, metric_palette, config
from utils.data_loader import fetch_data
from utils.preprocess import select_datasets
from app_pages import home, gene, allele, help, citations, custom_gene, custom_allele
from utils.auth import check_password

if not check_password():
    st.stop()

st.set_page_config(page_title="MWT Dashboard", layout="wide")
st.title("MWT Data Dashboard - Rankin Lab @ UBC")

# Select a page
page = st.sidebar.radio("Select a page", pages)

# Load and filter data
data = fetch_data()
if page not in [pages[5], pages[6]]: # If page is not "Help" or "Citations"
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
    help.render()
elif page == pages[6]:
    citations.render()


