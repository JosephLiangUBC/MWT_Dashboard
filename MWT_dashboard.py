import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import sqlite3
import io


# Fail gracefully option

@st.cache_data
def read(table):
    result = pd.read_sql_query(f"SELECT * FROM {table}", conn)
    return result


filepicker = st.file_uploader("Choose a .db file", accept_multiple_files=False)
filename = filepicker.name
# st.write(f"file picked is {filename}")

# Read data from SQLite database

conn = sqlite3.connect('/Users/Joseph/Desktop/NRSC510B/mwt_data.db')

tap_output = read('tap_response_data')
tap_tstat_allele = read('tstat_gene_data')
# allele_metric_data = read('allele_phenotype_data')
gene_profile_data = read('gene_profile_data')
allele_profile_data = read('allele_profile_data')
gene_MSD = read('gene_MSD')
allele_MSD = read('allele_MSD')

# st.write(gene_MSD[gene_MSD['Screen']=='G-Protein_Screen'])
# tap_output = pd.read_sql_query("SELECT * FROM tap_response_data", conn)
# tap_baseline = pd.read_sql_query("SELECT * FROM tap_baseline_data", conn)
conn.close()
tap_output['Strain'] = tap_output['Gene'] + " (" + tap_output['Allele'] + ")"

metric_palette = ["k", "k", "k",
                  "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgrey",
                  "lightsteelblue", "lightsteelblue", "lightsteelblue",
                  "powderblue", "powderblue", "powderblue",
                  "cadetblue", "cadetblue", "cadetblue",
                  "thistle", "thistle", "thistle"]

# Streamlit Dashboard starts here
st.title('NRSC510B: Data Dashboard for MWT Data')

datasets = st.multiselect(
    label="Select Datasets",
    options=gene_MSD.Screen.unique(),
    default=gene_MSD.Screen.unique()[0],
    placeholder="make a selection",
    help="select and de-select datasets you want to analyze",
    key="datasetselection"
)

phenotype_list = []
for a in gene_MSD.columns[1:]:
    b = a.split("-", 1)[0]
    phenotype_list.append(b)

dropna_features = list(np.unique(phenotype_list))
dropna_features.remove('Spontaneous Recovery of Response Duration')
dropna_features.remove('Spontaneous Recovery of Response Probability')
dropna_features.remove('Spontaneous Recovery of Response Speed')
# st.write(dropna_features)

tap_output = tap_output[tap_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
# tap_tstat_allele = tap_tstat_allele[tap_tstat_allele['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
tap_tstat_allele = tap_tstat_allele[tap_tstat_allele['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(
    columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
gene_profile_data = gene_profile_data[gene_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_profile_data = allele_profile_data[allele_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"],
                                                                                                "N2")
gene_MSD = gene_MSD[gene_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_MSD = allele_MSD[allele_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

#
# st.write(datasets)
# st.write(allele_profile_data)
# st.write(tap_output)
# st.write(gene_MSD)
st.header('Data at a glance')

# if st.checkbox('Show MSD data'):
#     st.subheader('Raw MSD data')
#     st.write(gene_MSD)

# st.write(gene_MSD[''])
# if st.checkbox('Show raw tap data'):
#     st.subheader('Raw tap data')
#     st.write(tap_output)
# #
# if st.checkbox('Show raw baseline data'):
#     st.subheader('Raw baseline data')
#     st.write(tap_baseline)
#
# if st.checkbox('Show tstat data'):
#     st.subheader('tstat data')
#     st.write(tap_tstat_allele)

col1, col2 = st.columns([2, 3])

col1.subheader("For A Single Phenotype")

phenotype_option = col1.selectbox(
    'Select a phenotype',
    np.unique(phenotype_list), key="phenotypeselect")

# seaborn graph of phenotypic view (sample mean distance if possible) + st.pyplot()
sns.set_context('notebook')
colors = ["dimgray"] * len(gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])["Gene"])
colors[gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]).reset_index(drop=True)[
    gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == "N2"].index[0]] = "red"
figx = col1.slider('Figure width', 1, 20, 4, key="figx_1")
figy = col1.slider('Figure Height', 1, 70, 16, key="figy_1")
fig, ax = plt.subplots(figsize=(figx, figy))
# fig, ax = plt.subplots()
# ax = sns.pointplot(data = gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]),
#             x=f"{phenotype_option}-mean",
#             y="Gene-",
#             # errorbar=list(zip(gene_MSD[f"{phenotype_option}-ci95_lo"],gene_MSD[f"{phenotype_option}-ci95_hi"])),
#             palette=["dimgray"]).set_title(f"{phenotype_option}")
ax = plt.errorbar(x=gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])[f"{phenotype_option}-mean"],
                  y=gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])["Gene"],
                  xerr=gene_MSD[f"{phenotype_option}-ci95_hi"] - gene_MSD[f"{phenotype_option}-mean"],
                  fmt="none", marker="none", ecolor=colors, elinewidth=3)
ax = plt.scatter(x=gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])[f"{phenotype_option}-mean"],
                 y=gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])["Gene"],
                 marker='o', color=colors)

plt.xlabel('Sample Mean Distance')
plt.ylabel('Gene')
plt.title(f"{phenotype_option}")
col1.pyplot(fig)
col1.caption(
    f'Sample mean distance from wildtype for all strains for selected phenotype: {phenotype_option}. Error bars are 95% CI')
phenotype_plot = io.BytesIO()
plt.savefig(phenotype_plot, format='png', dpi=300, bbox_inches='tight')
# Insert download graph button
col1.download_button(label="Download Plot",
                     data=phenotype_plot,
                     file_name=f"{phenotype_option}_profile.png",
                     mime="image/png",
                     key='dnldphenotypeprofile')

# Insert download graph button

col2.subheader("Comprehensive Heatmap")
sns.set_context('notebook', font_scale=0.7)
figx_hm = col2.slider('Figure Width', 0, 30, 15, key="figx_hm")
figy_hm = col2.slider('Figure Height', 0, 70, 20, key="figy_hm")
fig, ax = plt.subplots(figsize=(figx_hm, figy_hm))
# fig, ax = plt.subplots()
# ax = sns.heatmap(glue)

ax = sns.heatmap(data=tap_tstat_allele.set_index('Gene').drop(index="N2"),
                 annot=False,
                 linewidth=0.2,
                 square=False,
                 cmap="vlag",
                 #                  cmap=sns.diverging_palette(55, 250, s=100, l=40,as_cmap=True),
                 center=0,
                 vmax=3,
                 vmin=-3,
                 # xticklabels='auto',
                 # yticklabels='auto',
                 cbar_kws={"shrink": .05, "pad": 0.01})
ax.set(xlabel="", ylabel="")
ax.set_facecolor('xkcd:black')
col2.pyplot(fig)
col2.caption('Comprehensive heatmap of entire dataset')

imgheatmap = io.BytesIO()
plt.savefig(imgheatmap, format='png', dpi=300, bbox_inches='tight')
col2.download_button(label="Download Plot",
                     data=imgheatmap,
                     file_name="Heatmap.png",
                     mime="image/png",
                     key='dnldheatmap')
# Insert download graph button
# Insert download csv

st.header('Gene-specific Data')
gene_option = st.selectbox(
    'Select a gene',
    (tap_output['Gene'].unique()), key="geneselect")

tap_output_gene = tap_output[tap_output['Gene'] == gene_option]
# st.write(tap_output_allele)
# st.write(tap_output_allele['Date'].unique())
gene_tap_data = tap_output[tap_output['Date'].isin(tap_output_gene['Date'].unique())]
gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2', gene_option])]
gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)
# st.write(gene_tap_data_plot)


col3, col4 = st.columns([1, 1])
col3.subheader('phenotypic profile')

# seaborn plot
sns.set_context('notebook', font_scale=1)
fig, ax = plt.subplots(figsize=(5, 5))
ax = sns.barplot(x="Metric",  # <- Here we use seaborn as our graphing package.
                 y="T_score", orient='v',
                 data=gene_profile_data[gene_profile_data.Gene == f"{gene_option}"],
                 palette=metric_palette).set_title(f"{gene_option}")
plt.xticks(rotation=90)
plt.ylabel("Normalized T-Score")  # <- X-axis title
plt.ylim(-3, 3)

col3.pyplot(fig)

# download graph button
gene_profile_plot = io.BytesIO()
plt.savefig(gene_profile_plot, format='png', dpi=300, bbox_inches='tight')

col3.caption(f'Phenotypic profile of {gene_option}')
col3.download_button(label="Download Plot",
                     data=gene_profile_plot,
                     file_name=f"{gene_option}_profile.png",
                     mime="image/png",
                     key='dnldgeneprofile')

gene_phenotype_option = col4.selectbox(
    'Select a phenotype',
    np.unique(phenotype_list),
    key='gene_phenotype_select')

col4.subheader('Rank in phenotype')

# seaborn graph of phenotypic view (sample mean distance) + st.pyplot
sns.set_context('notebook')
gene_colors = ["dimgray"] * len(gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"])
gene_colors[gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
    gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == "N2"].index[
    0]] = "red"
gene_colors[gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
    gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == gene_option].index[
    0]] = "magenta"
fig, ax = plt.subplots(figsize=(figx, figy))
# ax = sns.pointplot(data = gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]),
#             x=f"{phenotype_option}-mean",
#             y="Gene-",
#             # errorbar=list(zip(gene_MSD[f"{phenotype_option}-ci95_lo"],gene_MSD[f"{phenotype_option}-ci95_hi"])),
#             palette=["dimgray"]).set_title(f"{phenotype_option}")
ax = plt.errorbar(x=gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-mean"],
                  y=gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"],
                  xerr=gene_MSD[f"{gene_phenotype_option}-ci95_hi"] - gene_MSD[f"{gene_phenotype_option}-mean"],
                  fmt="none", marker="none", ecolor=gene_colors, elinewidth=3)
ax = plt.scatter(x=gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-mean"],
                 y=gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"],
                 marker='o', color=gene_colors)

plt.xlabel('Sample Mean Distance')
plt.ylabel('Gene')
plt.title(f"{gene_phenotype_option}")

col4.pyplot(fig)

gene_phenotype_plot = io.BytesIO()
plt.savefig(gene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
col4.caption(
    f'Sample mean distance from wildtype for selected phenotype: {gene_phenotype_option}. Error bars are 95% CI.')
# Insert download graph button
col4.download_button(label="Download Plot",
                     data=gene_phenotype_plot,
                     file_name=f"{gene_option}_{gene_phenotype_option}_profile.png",
                     mime="image/png",
                     key='dnldgenephenotypeprofile')

st.subheader('Habituation Curves')
tab1, tab2, tab3 = st.tabs(["Habituation of Response Probability",
                            "Habituation of Response Duration",
                            "Habituation of Response Speed"])

with tab1:
    #  Habituation of Response Probability Plot
    st.subheader("Habituation of Response Probability")
    fig, ax = plt.subplots(figsize=(12, 10))
    # seaborn plot
    plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
    ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                       y="prob",
                       data=gene_tap_data_plot,
                       hue='Gene',  # <- Here we use the extra column from step 6 to separate by group
                       errorbar='se')  # <- Confidence interval. 95 = standard error
    plt.xlabel("Taps")  # <- X-axis title
    plt.ylabel("Probability")  # <- Y-Axis title
    plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
    plt.ylim(0, 1)
    ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

    tab1.pyplot(fig)
    tab1.caption(f'Habituation of Response Probability: {gene_option}')

    # download graph button
    img1 = io.BytesIO()
    plt.savefig(img1, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img1,
                       file_name=f"Probability of Tap Habituation {gene_option}.png",
                       mime="image/png",
                       key='dnldbtn1')

with tab2:
    #  Habituation of Response Duration Plot
    st.subheader("Habituation of Response Duration")
    fig, ax = plt.subplots(figsize=(12, 10))
    # seaborn plot
    ax = sns.pointplot(x="taps",
                       y="dura",
                       data=gene_tap_data_plot,
                       hue='Gene',
                       # palette=pal,
                       errorbar='se')
    plt.xlabel("Taps", fontsize='12')
    plt.ylabel("Duration", fontsize='12')
    plt.title("Habituation of Response Duration", fontsize='16')
    plt.ylim(0, None)
    ax.legend(loc='upper right', fontsize='12')
    tab2.pyplot(fig)
    tab2.caption(f'Habituation of Response Duration: {gene_option}')

    # download graph button
    img2 = io.BytesIO()
    plt.savefig(img2, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img2,
                       file_name=f"Duration of Tap Habituation {gene_option}.png",
                       mime="image/png",
                       key='dnldbtn2')
# Seaborn Graph of Duration Habituation curve
# Insert download graph button

with tab3:
    #  Habituation of Response Speed Plot
    st.subheader("Habituation of Response Speed")
    fig, ax = plt.subplots(figsize=(12, 10))
    # seaborn plot
    ax = sns.pointplot(x="taps",
                       y="speed",
                       data=gene_tap_data_plot,
                       hue='Gene',
                       errorbar='se')
    plt.xlabel("Taps", fontsize='12')
    plt.ylabel("Speed", fontsize='12')
    plt.title("Habituation of Response Speed", fontsize='16')
    plt.ylim(0, None)
    ax.legend(loc='upper right', fontsize='12')
    tab3.pyplot(fig)
    tab3.caption(f'Habituation of Response Speed: {gene_option}')

    # download graph button
    img3 = io.BytesIO()
    plt.savefig(img3, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img3,
                       file_name=f"Speed of Tap Habituation {gene_option}.png",
                       mime="image/png",
                       key='dnldbtn3')
# seaborn graph of Speed Habituation Curve
# Insert download graph button


st.header('Allele-specific Data')
allele_option = st.selectbox(
    'Select a allele',
    (tap_output['dataset'].unique()))

tap_output_allele = tap_output[tap_output['dataset'] == allele_option]
# st.write(tap_output_allele)
# st.write(tap_output_allele['Date'].unique())
allele_tap_data = tap_output[tap_output['Date'].isin(tap_output_allele['Date'].unique())]
allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2', allele_option])]
allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)
# st.write(allele_tap_data_plot)

col5, col6 = st.columns([1, 1])
col5.subheader('phenotypic profile')

# seaborn plot
sns.set_context('notebook', font_scale=1)
fig, ax = plt.subplots(figsize=(5, 5))
ax = sns.barplot(x="Metric",  # <- Here we use seaborn as our graphing package.
                 y="T_score", orient='v',
                 data=allele_profile_data[allele_profile_data.dataset == f"{allele_option}"],
                 palette=metric_palette).set_title(f"{allele_option}")
plt.xticks(rotation=90)

plt.ylabel("Normalized T-Score")  # <- X-axis title
plt.ylim(-3, 3)

col5.pyplot(fig)

# download graph button
allele_profile_plot = io.BytesIO()
plt.savefig(allele_profile_plot, format='png', dpi=300, bbox_inches='tight')

col5.caption(f'Phenotypic profile of gene-allele {allele_option}')
col5.download_button(label="Download Plot",
                     data=allele_profile_plot,
                     file_name=f"{allele_option}_profile.png",
                     mime="image/png",
                     key='dnldalleleprofile')

allele_phenotype_option = col6.selectbox(
    'Select a phenotype',
    np.unique(phenotype_list), key='allele_phenotype_select')
# seaborn graph of phenotypic view (sample mean distance) + st.pyplot

sns.set_context('notebook')
allele_colors = ["dimgray"] * len(allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])["dataset"])
allele_colors[allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"]).reset_index(drop=True)[
    allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"]).reset_index(drop=True)["dataset"] == "N2"].index[
    0]] = "red"
allele_colors[allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"]).reset_index(drop=True)[
    allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"]).reset_index(drop=True)[
        "dataset"] == allele_option].index[
    0]] = "magenta"

fig, ax = plt.subplots(figsize=(4, 16))
# ax = sns.pointplot(data = gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]),
#             x=f"{phenotype_option}-mean",
#             y="Gene-",
#             # errorbar=list(zip(gene_MSD[f"{phenotype_option}-ci95_lo"],gene_MSD[f"{phenotype_option}-ci95_hi"])),
#             palette=["dimgray"]).set_title(f"{phenotype_option}")
ax = plt.errorbar(x=allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])[f"{allele_phenotype_option}-mean"],
                  y=allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])["dataset"],
                  xerr=allele_MSD[f"{allele_phenotype_option}-ci95_hi"] - allele_MSD[f"{allele_phenotype_option}-mean"],
                  fmt="none", marker="none", ecolor=allele_colors, elinewidth=3)
ax = plt.scatter(x=allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])[f"{allele_phenotype_option}-mean"],
                 y=allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])["dataset"],
                 marker='o', color=allele_colors)

plt.xlabel('Sample Mean Distance')
plt.ylabel('Gene_Allele')
plt.title(f"{allele_phenotype_option}")

col6.pyplot(fig)

allele_phenotype_plot = io.BytesIO()
plt.savefig(allele_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
col6.caption(
    f'Sample mean distance from wildtype for selected phenotype: {allele_phenotype_option}. Error bars are 95% CI.')
# Insert download graph button
col6.download_button(label="Download Plot",
                     data=allele_phenotype_plot,
                     file_name=f"{allele_option}_{allele_phenotype_option}_profile.png",
                     mime="image/png",
                     key='dnldallelephenotypeprofile')

# Insert download graph button


st.subheader('Habituation Curves')
tab4, tab5, tab6 = st.tabs(["Habituation of Response Probability",
                            "Habituation of Response Duration",
                            "Habituation of Response Speed"])

with tab4:
    #  Habituation of Response Probability Plot
    st.subheader("Habituation of Response Probability")
    fig, ax = plt.subplots(figsize=(12, 10), linewidth=2.5)
    # seaborn plot
    plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
    ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                       y="prob",
                       data=allele_tap_data_plot,
                       hue='dataset',  # <- Here we use the extra column from step 6 to separate by group
                       errorbar='se')  # <- Confidence interval. 95 = standard error
    plt.xlabel("Taps")  # <- X-axis title
    plt.ylabel("Probability")  # <- Y-Axis title
    plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
    plt.ylim(0, 1)
    ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

    tab4.pyplot(fig)
    tab4.caption(f'Habituation of Response Probability: {allele_option}')

    # download graph button
    img4 = io.BytesIO()
    plt.savefig(img4, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img4,
                       file_name=f"Probability of Tap Habituation {allele_option}.png",
                       mime="image/png",
                       key='dnldbtn4')

with tab5:
    #  Habituation of Response Duration Plot
    st.subheader("Habituation of Response Duration")
    fig, ax = plt.subplots(figsize=(12, 10))
    # seaborn plot
    ax = sns.pointplot(x="taps",
                       y="dura",
                       data=allele_tap_data_plot,
                       hue='dataset',
                       # palette=pal,
                       errorbar='se')
    plt.xlabel("Taps", fontsize='12')
    plt.ylabel("Duration", fontsize='12')
    plt.title("Habituation of Response Duration", fontsize='16')
    plt.ylim(0, None)
    ax.legend(loc='upper right', fontsize='12')
    tab5.pyplot(fig)
    tab5.caption(f'Habituation of Response Duration: {allele_option}')

    # download graph button
    img5 = io.BytesIO()
    plt.savefig(img5, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img5,
                       file_name=f"Duration of Tap Habituation {allele_option}.png",
                       mime="image/png",
                       key='dnldbtn5')

with tab6:
    #  Habituation of Response Speed Plot
    st.subheader("Habituation of Response Speed")
    fig, ax = plt.subplots(figsize=(12, 10))
    # seaborn plot
    ax = sns.pointplot(x="taps",
                       y="speed",
                       data=allele_tap_data_plot,
                       hue='dataset',
                       errorbar='se')
    plt.xlabel("Taps", fontsize='12')
    plt.ylabel("Speed", fontsize='12')
    plt.title("Habituation of Response Speed", fontsize='16')
    plt.ylim(0, None)
    ax.legend(loc='upper right', fontsize='12')
    tab6.pyplot(fig)
    tab6.caption(f'Habituation of Response Speed: {allele_option}')

    # download graph button
    img6 = io.BytesIO()
    plt.savefig(img6, format='png', dpi=300, bbox_inches='tight')
    st.download_button(label="Download Plot",
                       data=img6,
                       file_name=f"Speed of Tap Habituation {allele_option}.png",
                       mime="image/png",
                       key='dnldbtn6')
