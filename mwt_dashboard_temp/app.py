import streamlit as st
from config import pages, metric_palette, config
from utils.data_loader import fetch_data
from utils.preprocess import select_datasets
from pages import home, gene, allele  # import other pages as needed
from utils.auth import check_password

# if not check_password():
#     st.stop()

st.set_page_config(page_title="MWT Dashboard", layout="wide")
st.title("MWT Data Dashboard - Rankin Lab @ UBC")

# Select a page
page = st.sidebar.radio("Select a page", pages)

# Load and filter data
data = fetch_data()
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
# Add more pages as needed

