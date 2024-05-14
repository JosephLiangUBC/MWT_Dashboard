import random
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import sqlite3
import io

# Added configuration to wide for desktop screen (do not remove)
# it has to be the first command in the file for it to run
# Streamlit command starts here
st.set_page_config(layout="wide")

# convert dataframe to csv for download
@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")

# Fail gracefully option
def read(table):
    result = pd.read_sql_query(f"SELECT * FROM {table}", conn)
    return result

conn = sqlite3.connect('/Users/Joseph/Desktop/NRSC510B/mwt_data.db')
# conn = sqlite3.connect('/Users/lavanya/Downloads/MWT_Dashboard-main/Test/mwt_data.db')

# Read data from SQLite database
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

# Defining the color palette for the plots
metric_palette = ["k", "k", "k",
                  "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgrey",
                  "lightsteelblue", "lightsteelblue", "lightsteelblue",
                  "powderblue", "powderblue", "powderblue",
                  "cadetblue", "cadetblue", "cadetblue",
                  "thistle", "thistle", "thistle"]

# Streamlit Dashboard title
st.title('NRSC510B: Data Dashboard for MWT Data')

# Select dataset option
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

# filter data for selected dataset
tap_output = tap_output[tap_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
# tap_tstat_allele = tap_tstat_allele[tap_tstat_allele['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
tap_tstat_allele = tap_tstat_allele[tap_tstat_allele['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(
    columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
gene_profile_data = gene_profile_data[gene_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_profile_data = allele_profile_data[allele_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"],
                                                                                                "N2")
gene_MSD = gene_MSD[gene_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_MSD = allele_MSD[allele_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

# creating tabs for dashboard
tabs_font_css = """
<style>
    .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
    font-size:24px;
    }
</style>
"""
st.write(tabs_font_css, unsafe_allow_html=True)
data_tab, gene_tab, allele_tab , custom_select_tab,clustering_tab= st.tabs(["Data at a Glance", "Gene-specific Data", "Allele-specific Data",  "Custom Selection","Clustering"])

# Visualisations for data tab
with data_tab:
    st.header('Data at a glance')

#
# st.write(datasets)
# st.write(allele_profile_data)
# st.write(tap_output)
# st.write(gene_MSD)

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
    fig, ax = plt.subplots(figsize=(4,16))
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
    
    phenotype_plot = io.BytesIO()
    plt.savefig(phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col1.image(phenotype_plot, width=None, caption=(f'Sample mean distance from wildtype for all strains for selected phenotype: {phenotype_option}. Error bars are 95% CI')) ## added
        
    #combine data and rename columns :
    data_dat=pd.concat([gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])["Gene"],
                        gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])[f"{phenotype_option}-mean"],
                        gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])[f"{phenotype_option}-ci95_lo"],
                        gene_MSD.sort_values(by=[f"{phenotype_option}-mean"])[f"{phenotype_option}-ci95_hi"]], 
                        axis=1)
    data_dat.columns=["Gene", f"{phenotype_option}", f"{phenotype_option}-lower" ,f"{phenotype_option}-upper"]

    # Insert download graph button
    col1.download_button(label="Download Plot",
                        data=phenotype_plot,
                        file_name=f"{phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldphenotypeprofile')
    col1.download_button(label="Download csv",
                            data=convert_df(data_dat),
                            file_name=f"Data Glance Sample mean distance {phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldphenotypeprofilecsv')

    # Insert download graph button

    col2.subheader("Comprehensive Heatmap")
    sns.set_context('notebook', font_scale=0.7)
    fig, ax = plt.subplots(figsize=(15, 20))
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

    imgheatmap = io.BytesIO()
    plt.savefig(imgheatmap, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col2.image(imgheatmap,caption='Comprehensive heatmap of entire dataset', width=None)
    col2.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmap')
    col2.download_button(label="Download csv",
                            data=convert_df(tap_tstat_allele.set_index('Gene').drop(index="N2")),
                            file_name=f"Data Glance Heatmap.csv",
                            mime="text/csv",
                            key='dnldheatmapcsv')
    # Insert download graph button
    # Insert download csv
with gene_tab:
    st.header('Gene-specific Data')
    gene_option = st.selectbox(
        'Select a gene',
        (tap_output['Gene'].unique()), key="geneselect")
    url_data= pd.read_csv("WB_id.csv")
    if gene_option:
        gene_id=url_data[url_data['Gene']==gene_option]['Identifier'].values[0]
        glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
    st.markdown(f'<p style="font-size:20px">For more information on <a href="{glink}">{gene_option}</a></p>', unsafe_allow_html=True)

    tap_output_gene = tap_output[tap_output['Gene'] == gene_option]
    # st.write(tap_output_allele)
    # st.write(tap_output_allele['Date'].unique())
    gene_tap_data = tap_output[tap_output['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2', gene_option])]
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)
    # st.write(gene_tap_data_plot)

 # added extra column to show abituation curves in landscape 
    col3, col4, col7 = st.columns([1, 1, 1])
    col3.subheader('Phenotypic profile')

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

    # download graph button
    gene_profile_plot = io.BytesIO()
    plt.savefig(gene_profile_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col3.image(gene_profile_plot,width=None,caption=(f'Phenotypic profile of {gene_option}.'))
    col3.download_button(label="Download Plot",
                        data=gene_profile_plot,
                        file_name=f"{gene_option}_profile.png",
                        mime="image/png",
                        key='dnldgeneprofile')
    col3.download_button(label="Download csv",
                        data=convert_df(gene_profile_data[gene_profile_data.Gene == f"{gene_option}"]),
                        file_name=f"Phenotypic profile of {gene_option}.csv",
                        mime="text/csv",
                        key='dnldgeneprofilecsv')
    
    col4.subheader('Rank in phenotype')
    gene_phenotype_option = col4.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list),
        key='gene_phenotype_select')


    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    sns.set_context('notebook')
    gene_colors = ["dimgray"] * len(gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"])
    gene_colors[gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == "N2"].index[
        0]] = "red"
    gene_colors[gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == gene_option].index[
        0]] = "magenta"
    fig, ax = plt.subplots(figsize=(4, 16))
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
    plt.yticks(fontsize=7) # added to see the axis labels better

    plt.xlabel('Sample Mean Distance')
    plt.ylabel('Gene')
    plt.title(f"{gene_phenotype_option}")

    gene_phenotype_plot = io.BytesIO()
    plt.savefig(gene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col4.image(gene_phenotype_plot, width=None,caption=f'Sample mean distance from wildtype for selected phenotype: {gene_phenotype_option}. Error bars are 95% CI.')
    
    #combine data and rename columns :
    gene_dat=pd.concat( [gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-mean"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-ci95_lo"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-ci95_hi"]],
                   axis=1)
    gene_dat.columns=["Gene", f"{gene_phenotype_option}", f"{gene_phenotype_option}-lower" ,f"{gene_phenotype_option}-upper"]
    
    # Insert download graph button
    col4.download_button(label="Download Plot",
                        data=gene_phenotype_plot,
                        file_name=f"{gene_option}_{gene_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldgenephenotypeprofile')
    col4.download_button(label="Download csv",
                            data=convert_df(gene_dat),
                            file_name=f"Gene-specific Data Sample mean distance {gene_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldgenephenotypeprofilecsv')

    # Move tabs inside col 7 for better viewing
    col7.subheader('Habituation Curves')
    with col7:
        tab1, tab2, tab3 = st.tabs(["Habituation of Response Probability",
                                "Habituation of Response Duration",
                                "Habituation of Response Speed"])
        with tab1:
            #  Habituation of Response Probability Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
            ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                            y="prob",
                            data=gene_tap_data_plot,
                            hue='Gene',  # <- Here we use the extra column from step 6 to separate by group
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
                            errorbar='se')  # <- Confidence interval. 95 = standard error
            plt.xlabel("Taps")  # <- X-axis title
            plt.ylabel("Probability")  # <- Y-Axis title
            plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
            plt.ylim(0, 1)
            ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

            # download graph button
            img1 = io.BytesIO()
            plt.savefig(img1, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab1.image(img1, width=None,caption=(f'Habituation of Response Probability: {gene_option}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img1,
                            file_name=f"Probability of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn1')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_option}.csv",
                            mime="text/csv",
                            key='dnldbtn7')

        with tab2:
            #  Habituation of Response Duration Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="dura",
                            data=gene_tap_data_plot,
                            hue='Gene',
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Duration", fontsize='12')
            plt.title("Habituation of Response Duration", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
            
            # download graph button
            img2 = io.BytesIO()
            plt.savefig(img2, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab2.image(img2, width=None,caption=(f'Habituation of Response Duration: {gene_option}'))
            # Insert download plot and download csv button

            st.download_button(label="Download Plot",
                            data=img2,
                            file_name=f"Duration of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn2')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_option}.csv",
                            mime="text/csv",
                            key='dnldbtn8')
        # Seaborn Graph of Duration Habituation curve

        with tab3:
            #  Habituation of Response Speed Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="speed",
                            data=gene_tap_data_plot,
                            hue='Gene',
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Speed", fontsize='12')
            plt.title("Habituation of Response Speed", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
        
            img3 = io.BytesIO()
            plt.savefig(img3, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab3.image(img3, width=None,caption=(f'Habituation of Response Speed: {gene_option}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img3,
                            file_name=f"Speed of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn3')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_option}.csv",
                            mime="text/csv",
                            key='dnldbtn9')
        # seaborn graph of Speed Habituation Curve
        # Insert download graph button

with allele_tab:
    st.header('Allele-specific Data')
    # select allele 
    allele_option = st.selectbox(
        'Select a allele',
        (tap_output['dataset'].unique()))
    #     wlink=f'https://wormbase.org/species/c_elegans/gene/{allele_id}'
    # st.markdown(f'For more allele information on [{allele_option}](%s)' % wlink)
    split=allele_option.split('_')
    url_data= pd.read_csv("WB_id.csv")
    if allele_option:
        gene_id=url_data[url_data['Gene']==split[0]]['Identifier'].values[0]
        glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
    st.markdown(f'<p style="font-size:20px">For more gene information on <a href="{glink}">{split[0]}</a></p>', unsafe_allow_html=True)

    tap_output_allele = tap_output[tap_output['dataset'] == allele_option]
    # st.write(tap_output_allele)
    # st.write(tap_output_allele['Date'].unique())
    allele_tap_data = tap_output[tap_output['Date'].isin(tap_output_allele['Date'].unique())]
    allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2', allele_option])]
    allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)
    # st.write(allele_tap_data_plot)

    col5, col6, col8 = st.columns([1, 1, 1])
    col5.subheader('Phenotypic profile')

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

    # download graph button
    allele_profile_plot = io.BytesIO()
    plt.savefig(allele_profile_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col5.image(allele_profile_plot,width=None,caption=(f'Phenotypic profile of gene-allele {allele_option}'))
    col5.download_button(label="Download Plot",
                        data=allele_profile_plot,
                        file_name=f"{allele_option}_profile.png",
                        mime="image/png",
                        key='dnldalleleprofile')
    col5.download_button(label="Download csv",
                        data=convert_df(allele_profile_data[allele_profile_data.dataset == f"{allele_option}"]),
                        file_name=f"Phenotypic profile of gene-allele {allele_option}.csv",
                        mime="text/csv",
                        key='dnldalleleprofilecsv')
    col6.subheader('Rank in phenotype')

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
    plt.yticks(fontsize=5)

    plt.xlabel('Sample Mean Distance')
    plt.ylabel('Gene_Allele')
    plt.title(f"{allele_phenotype_option}")

    allele_phenotype_plot = io.BytesIO()
    plt.savefig(allele_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col6.image(allele_phenotype_plot,width=None,caption=(f'Sample mean distance from wildtype for selected phenotype: {allele_phenotype_option}. Error bars are 95% CI.'))

    #combine data and rename columns :
    dat=pd.concat( [allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])["dataset"],
                   allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])[f"{allele_phenotype_option}-mean"],
                   allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])[f"{allele_phenotype_option}-ci95_lo"],
                   allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])[f"{allele_phenotype_option}-ci95_hi"]],
                   axis=1)
    dat.columns=["gene-allele", f"{allele_phenotype_option}", f"{allele_phenotype_option}-lower" ,f"{allele_phenotype_option}-upper"]
    
    # Insert download graph button
    col6.download_button(label="Download Plot",
                        data=allele_phenotype_plot,
                        file_name=f"{allele_option}_{allele_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldallelephenotypeprofile')
    col6.download_button(label="Download csv",
                        data=convert_df(dat),
                        file_name=f"Allele-specific Data Sample mean distance {allele_phenotype_option}.csv",
                        mime="text/csv",
                        key='dnldallelephenotypeprofilecsv')
    # Insert download graph button


    col8.subheader('Habituation Curves')
    with col8:
        tab4, tab5, tab6 = st.tabs(["Habituation of Response Probability",
                                    "Habituation of Response Duration",
                                    "Habituation of Response Speed"])

        with tab4:
            #  Habituation of Response Probability Plot
            fig, ax = plt.subplots(figsize=(12, 10), linewidth=2.5)
            # seaborn plot
            plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
            ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                            y="prob",
                            data=allele_tap_data_plot,
                            hue='dataset',  # <- Here we use the extra column from step 6 to separate by group
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
                            errorbar='se')  # <- Confidence interval. 95 = standard error
            plt.xlabel("Taps")  # <- X-axis title
            plt.ylabel("Probability")  # <- Y-Axis title
            plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
            plt.ylim(0, 1)
            ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

            img4 = io.BytesIO()
            plt.savefig(img4, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab4.image(img4, width=None,caption=(f'Habituation of Response Probability: {allele_option}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img4,
                            file_name=f"Probability of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn4')
            st.download_button(label="Download csv",
                                data=convert_df(allele_tap_data_plot),
                                file_name=f"Allele-specific Data {allele_option}.csv",
                                mime="text/csv",
                                key='dnldbtn10')

        with tab5:
            #  Habituation of Response Duration Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="dura",
                            data=allele_tap_data_plot,
                            hue='dataset',
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Duration", fontsize='12')
            plt.title("Habituation of Response Duration", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
            img5 = io.BytesIO()
            plt.savefig(img5, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab5.image(img5, width=None,caption=(f'Habituation of Response Duration: {allele_option}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img5,
                            file_name=f"Duration of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn5')
            st.download_button(label="Download csv",
                                data=convert_df(allele_tap_data_plot),
                                file_name=f"Allele-specific Data {allele_option}.csv",
                                mime="text/csv",
                                key='dnldbtn11')

        with tab6:
            #  Habituation of Response Speed Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="speed",
                            data=allele_tap_data_plot,
                            hue='dataset',
                            palette=["steelblue" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Speed", fontsize='12')
            plt.title("Habituation of Response Speed", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
            
            img6 = io.BytesIO()
            plt.savefig(img6, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab6.image(img6, width=None,caption=(f'Habituation of Response Speed: {allele_option}'))        
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img6,
                            file_name=f"Speed of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn6')
            st.download_button(label="Download csv",
                                data=convert_df(allele_tap_data_plot),
                                file_name=f"Allele-specific Data {allele_option}.csv",
                                mime="text/csv",
                                key='dnldbtn12')

with custom_select_tab:
   # multiple selection option for genes
    gene_multiple = st.multiselect(
    label="Select Genes",
    options=tap_output['Gene'].unique(),
    default=tap_output['Gene'].unique()[0],
    placeholder="make a selection",
    help="select and de-select genes you want to analyze",
    key="geneselection")
  
    #fiilter data for particular genes
    tap_output_gene = tap_output[tap_output['Gene'].isin(gene_multiple)]
    # st.write(tap_output_allele)
    # st.write(tap_output_allele['Date'].unique())
    gene_tap_data = tap_output[tap_output['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2']+ gene_multiple)]
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)
    
    #add columns for msd, habituation plots and heatmap plots
    col9, col10, col11= st.columns([1,1,1])
    col9.subheader('Rank in phenotype')
    multigene_phenotype_option = col9.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list),
        key='multigene_phenotype_select')
    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    sns.set_context('notebook')
    gene_colors = ["dimgray"] * len(gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"])
    gene_colors[gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == "N2"].index[
        0]] = "red"
    for g in gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"].isin(gene_multiple)].index:
        gene_colors[g] = "magenta"
    fig, ax = plt.subplots(figsize=(4, 16))
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
    plt.yticks(fontsize=7) # added to see the axis labels better

    plt.xlabel('Sample Mean Distance')
    plt.ylabel('Genes')
    plt.title(f"{gene_phenotype_option}")

    multigene_phenotype_plot = io.BytesIO()
    plt.savefig(multigene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col9.image(multigene_phenotype_plot, width=None,caption=f'Sample mean distance from wildtype for selected phenotype: {gene_phenotype_option} and selected genes :{gene_multiple}. Error bars are 95% CI.')
    
    #combine data and rename columns :
    multigene_dat=pd.concat( [gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])["Gene"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-mean"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-ci95_lo"],
                   gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])[f"{gene_phenotype_option}-ci95_hi"]],
                   axis=1)
    gene_dat.columns=["Gene", f"{gene_phenotype_option}", f"{gene_phenotype_option}-lower" ,f"{gene_phenotype_option}-upper"]
    
    # Insert download graph button
    col9.download_button(label="Download Plot",
                        data=multigene_phenotype_plot,
                        file_name=f"multi_gene_{gene_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldmultigenephenotypeprofile')
    col9.download_button(label="Download csv",
                            data=convert_df(multigene_dat),
                            file_name=f"Gene-specific Data Sample mean distance {gene_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldmultigenephenotypeprofilecsv')
    
    col10.subheader('Habituation Curves')
    genes = gene_tap_data_plot['Gene'].unique()

    # Create a list of unique colors for the genes
    colors = []
    while len(colors) < len(genes):
        color = sns.color_palette()[random.randint(1, len(sns.color_palette())-1)]
        if color not in colors:
            colors.append(color)

    # Create a palette with 'teelblue' for 'N2' and the unique colors for the other genes
    new_palette = ["steelblue" if gene == "N2" else color for gene, color in zip(genes, colors)]

    with col10:
        tab7, tab8, tab9 = st.tabs(["Habituation of Response Probability",
                                "Habituation of Response Duration",
                                "Habituation of Response Speed"])
        with tab7:
            #  Habituation of Response Probability Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
            ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                            y="prob",
                            data=gene_tap_data_plot,
                            hue='Gene',  # <- Here we use the extra column from step 6 to separate by group
                            palette=new_palette,
                            errorbar='se')  # <- Confidence interval. 95 = standard error
            plt.xlabel("Taps")  # <- x-axis title
            plt.ylabel("Probability")  # <- y-axis title
            plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
            plt.ylim(0, 1)
            ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

            # download graph button
            img7 = io.BytesIO()
            plt.savefig(img7, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab7.image(img7, width=None,caption=(f'Habituation of Response Probability: {gene_multiple}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img7,
                            file_name=f"Probability of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn13')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn14')

        with tab8:
            #  Habituation of Response Duration Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="dura",
                            data=gene_tap_data_plot,
                            hue='Gene',
                            palette=new_palette, # N2 to be blue consistently
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Duration", fontsize='12')
            plt.title("Habituation of Response Duration", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
            
            # download graph button
            img8 = io.BytesIO()
            plt.savefig(img8, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab8.image(img8, width=None,caption=(f'Habituation of Response Duration: {gene_multiple}'))
            # Insert download plot and download csv button

            st.download_button(label="Download Plot",
                            data=img8,
                            file_name=f"Duration of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn15')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn16')
        # Seaborn Graph of Duration Habituation curve

        with tab9:
            #  Habituation of Response Speed Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="speed",
                            data=gene_tap_data_plot,
                            hue='Gene',
                            palette=new_palette, # N2 to be blue consistently
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Speed", fontsize='12')
            plt.title("Habituation of Response Speed", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
        
            img9 = io.BytesIO()
            plt.savefig(img9, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab9.image(img9, width=None,caption=(f'Habituation of Response Speed: {gene_multiple}'))
            # Insert download plot and download csv button
            st.download_button(label="Download Plot",
                            data=img9,
                            file_name=f"Speed of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn17')
            st.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn18')
        # seaborn graph of Speed Habituation Curve
        # Insert download graph button

    col11.subheader("Comprehensive Heatmap")
    sns.set_context('notebook', font_scale=0.7)
    fig, ax = plt.subplots(figsize=(15, 20))
    # fig, ax = plt.subplots()
    # ax = sns.heatmap(glue)
    # Filter the dataframe for the selected genes
    tap_tstat_allele_selected = tap_tstat_allele[tap_tstat_allele['Gene'].isin(gene_multiple)]

    ax = sns.heatmap(data=tap_tstat_allele_selected.set_index('Gene'),
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

    imgheatmap = io.BytesIO()
    plt.savefig(imgheatmap, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col11.image(imgheatmap,caption='Comprehensive heatmap of the dataset with selected genes', width=None)
    col11.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmapcustom')
    col11.download_button(label="Download csv",
                            data=convert_df(tap_tstat_allele_selected.set_index('Gene')),
                            file_name=f"Data Glance Heatmap.csv",
                            mime="text/csv",
                            key='dnldheatmapcsvcustom')

