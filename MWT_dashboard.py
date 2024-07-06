import itertools
import random
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import sqlite3
import io
import plotly.graph_objects as go

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

# list of connections to run the dashboard
conn_list=['/Users/Joseph/Desktop/NRSC510B/mwt_data.db',
           '/Users/lavanya/Downloads/MWT_Dashboard-main/Test/mwt_data_updated.db',
           '/Users/rankinlab/Desktop/MWT_Data_App/mwt_data.db']
for conn_path in conn_list:
    try:
        conn= sqlite3.connect(conn_path) #tries whichever connection is available
        break
    except:
        pass
        
# Read data from SQLite database
tap_output = read('tap_response_data')
baseline_output=read('tap_baseline_data') # for raw baseline data
tap_tstat_allele = read('tstat_gene_data')
tap_tstat_data = read('tstat_allele_data')
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
st.title('Data Dashboard for MWT Data')

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

phenotype_list.remove('Screen')

dropna_features = list(np.unique(phenotype_list))
dropna_features.remove('Spontaneous Recovery of Response Duration')
dropna_features.remove('Spontaneous Recovery of Response Probability')
dropna_features.remove('Spontaneous Recovery of Response Speed')
# st.write(dropna_features)

# filter data for selected dataset
tap_output = tap_output[tap_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
tap_tstat_allele = tap_tstat_allele[tap_tstat_allele['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(
    columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
tap_tstat_data = tap_tstat_data[tap_tstat_data['Screen'].isin(datasets)].dropna(subset=dropna_features).drop(
    columns=['Screen']).replace(["N2_N2", "N2_XJ1"], "N2")
gene_profile_data = gene_profile_data[gene_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_profile_data = allele_profile_data[allele_profile_data['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"],
                                                                                                "N2")
gene_MSD = gene_MSD[gene_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
allele_MSD = allele_MSD[allele_MSD['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")

# creating tabs for dashboard
pages = ["Data at a Glance", "Gene-specific Data", "Allele-specific Data",  "Custom Gene Selection","Custom Allele Selection","Clustering"]
page = st.sidebar.radio("Select a page", pages)

# Visualisations for data tab
if page ==pages[0]:
    st.header('Data at a glance')

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
    col1_1,col1_2=col1.columns(2)
    col1_1.download_button(label="Download Plot",
                        data=phenotype_plot,
                        file_name=f"{phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldphenotypeprofile')
    col1_2.download_button(label="Download csv",
                            data=convert_df(data_dat),
                            file_name=f"Data Glance Sample mean distance {phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldphenotypeprofilecsv')

    # Insert download graph button

    # Create a heatmap
    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_allele.set_index('Gene').drop(index="N2").values,
        x=tap_tstat_allele.set_index('Gene').drop(index="N2").columns,
        y=tap_tstat_allele.set_index('Gene').drop(index="N2").index,
        colorscale='RdBu',
        zmin=-3,
        zmax=3,
        colorbar=dict(
            len=0.95,
            thickness=10,
            tickvals=[-3, 0, 3],
            ticktext=['-3', '0', '3'],
            title="",
            titleside="right"
        )
    ))

    fig.update_layout(
        width=900,
        height=1200,
        margin=dict(l=50, r=50, t=50, b=50),
        xaxis_title="",
        yaxis_title=""
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col2.subheader("Comprehensive Heatmap of entire dataset")
    col2.plotly_chart(fig, use_container_width=True)

    # Add download buttons    
    col2.download_button(
        label="Download CSV",
        data=convert_df(tap_tstat_allele.set_index('Gene').drop(index="N2")),
        file_name="Data_Glance_Heatmap.csv",
        mime="text/csv",
        key='dnldheatmapcsv'
    )

if page == pages[1]:
    st.header('Gene-specific Data')
    gene_option = st.selectbox(
        'Select a gene',
        [gene for gene in tap_output['Gene'].unique()if gene != 'N2'], 
        key="geneselect")
    gene_id_data= pd.read_csv("WB_id.csv")

    if gene_option:
        try:
            gene_id=gene_id_data.loc[gene_id_data['Gene']==gene_option,'Identifier'].values[0]
            glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
            st.markdown(f'<p style="font-size:20px">For more gene information on <a href="{glink}">{gene_option}</a> (Source: GenomeAlliance)</p>', unsafe_allow_html=True)
        except:
            glink = 'https://www.alliancegenome.org'
            st.markdown(f"information for {gene_option} not available")
    

    tap_output_gene = tap_output[tap_output['Gene'] == gene_option]
    gene_tap_data = tap_output[tap_output['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2', gene_option])].dropna(subset=['taps'])
    # st.write(gene_tap_data_plot)
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)

 # added extra column to show habituation curves in landscape 
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
    col3_1,col3_2= col3.columns(2)
    col3_1.download_button(label="Download Plot",
                        data=gene_profile_plot,
                        file_name=f"{gene_option}_profile.png",
                        mime="image/png",
                        key='dnldgeneprofile')
    col3_2.download_button(label="Download csv",
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
    col4_1,col4_2= col4.columns(2)
    col4_1.download_button(label="Download Plot",
                        data=gene_phenotype_plot,
                        file_name=f"{gene_option}_{gene_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldgenephenotypeprofile')
    col4_2.download_button(label="Download csv",
                            data=convert_df(gene_dat),
                            file_name=f"Gene-specific Data Sample mean distance {gene_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldgenephenotypeprofilecsv')

    # Move tabs inside col 7 for better viewing
    col7.subheader('Habituation Curves of Response')
    with col7:
        tab1, tab2, tab3 = st.tabs(["Probability",
                                "Duration",
                                "Speed"])
        
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
            tab1_1,tab1_2= tab1.columns(2)

            tab1_1.download_button(label="Download Plot",
                            data=img1,
                            file_name=f"Probability of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn1')
            tab1_2.download_button(label="Download csv",
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
            tab2_1,tab2_2=tab2.columns(2)
            tab2_1.download_button(label="Download Plot",
                            data=img2,
                            file_name=f"Duration of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn2')
            tab2_2.download_button(label="Download csv",
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
            tab3_1,tab3_2=tab3.columns(2)
            tab3_1.download_button(label="Download Plot",
                            data=img3,
                            file_name=f"Speed of Tap Habituation {gene_option}.png",
                            mime="image/png",
                            key='dnldbtn3')
            tab3_2.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_option}.csv",
                            mime="text/csv",
                            key='dnldbtn9')
        # seaborn graph of Speed Habituation Curve
        # Insert download graph button
    st.download_button(label="Download raw baseline data",
                    data=convert_df(baseline_output),
                    file_name=f"raw_baseline_data.csv",
                    mime="text/csv",
                    key='dnldgenebaseoutcsv')

if page ==pages[2]:
    st.header('Allele-specific Data')
    # select allele 
    allele_option = st.selectbox(
        'Select a allele',
        [allele for allele in tap_output['dataset'].unique() if allele != 'N2'])
    #splititing gene allele to get gene and allele columns
    gene, allele = allele_option.split('_')
    gene_id_data= pd.read_csv("WB_id.csv")# data for gene id
    allele_id_data= pd.read_csv("Gene_Allele_WormBaseID.csv",names=['Identifier', 'Gene', 'Allele'])# data for allele id

    # check if allele and gene options match the preexisting list
    if allele_option:
        gene_id = gene_id_data.loc[gene_id_data['Gene'] == gene, 'Identifier'].values[0]
        allele_id = allele_id_data.loc[(allele_id_data['Gene'] == gene) & (allele_id_data['Allele'] == allele), 'Identifier']
        if allele_id.any():
            allele_id=allele_id.values[0]
        else: # if allele option doesnt match then send to default page
            allele_id= "default value"

        # links to the Alliance Genome and WormBase website 
        glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
        wlink=f'https://wormbase.org/species/c_elegans/variation/{allele_id}'
    # display links
    st.markdown(f'<p style="font-size:20px">For more gene information on <a href="{glink}">{gene}</a>(Source: GenomeAlliance) and allele information on <a href="{wlink}">{allele}</a>(Source: WormBase)</p>', unsafe_allow_html=True)

    tap_output_allele = tap_output[tap_output['dataset'] == allele_option]
    allele_tap_data = tap_output[tap_output['Date'].isin(tap_output_allele['Date'].unique())]
    allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2', allele_option])].dropna(subset=['taps'])
    allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)

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
    # display image 
    col5.image(allele_profile_plot,width=None,caption=(f'Phenotypic profile of gene-allele {allele_option}'))
    # download button for plots
    col5_1,col5_2=col5.columns(2)
    col5_1.download_button(label="Download Plot",
                        data=allele_profile_plot,
                        file_name=f"{allele_option}_profile.png",
                        mime="image/png",
                        key='dnldalleleprofile')
    
    # download button for data
    col5_2.download_button(label="Download csv",
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
    col6_1,col6_2=col6.columns(2)
    col6_1.download_button(label="Download Plot",
                        data=allele_phenotype_plot,
                        file_name=f"{allele_option}_{allele_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldallelephenotypeprofile')
    col6_2.download_button(label="Download csv",
                        data=convert_df(dat),
                        file_name=f"Allele-specific Data Sample mean distance {allele_phenotype_option}.csv",
                        mime="text/csv",
                        key='dnldallelephenotypeprofilecsv')
    # Insert download graph button


    col8.subheader('Habituation Curves of Response')
    with col8:
        tab4, tab5, tab6 = st.tabs(["Probability",
                                    "Duration",
                                    "Speed"])

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
            tab4_1,tab4_2=tab4.columns(2)
            tab4_1.download_button(label="Download Plot",
                            data=img4,
                            file_name=f"Probability of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn4')
            tab4_2.download_button(label="Download csv",
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
            tab5_1,tab5_2=tab5.columns(2)
            tab5_1.download_button(label="Download Plot",
                            data=img5,
                            file_name=f"Duration of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn5')
            tab5_2.download_button(label="Download csv",
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
            tab6_1,tab6_2=tab6.columns(2)
            tab6_1.download_button(label="Download Plot",
                            data=img6,
                            file_name=f"Speed of Tap Habituation {allele_option}.png",
                            mime="image/png",
                            key='dnldbtn6')
            tab6_2.download_button(label="Download csv",
                                data=convert_df(allele_tap_data_plot),
                                file_name=f"Allele-specific Data {allele_option}.csv",
                                mime="text/csv",
                                key='dnldbtn12')
    st.download_button(label="Download raw baseline data",
                    data=convert_df(baseline_output),
                    file_name=f"raw_baseline_data.csv",
                    mime="text/csv",
                    key='dnldallelebaseoutcsv')

if page ==pages[3]:
   # multiple selection option for genes
    st.header('Custom Gene Selection ')

    gene_multiple = st.multiselect(
        label="Select Genes",
        options=[gene for gene in tap_output['Gene'].unique() if gene != 'N2'],
        default=[gene for gene in tap_output['Gene'].unique() if gene != 'N2'][0],
        placeholder="make a selection",
        help="select and de-select genes you want to analyze",
        key="geneselection")
    gene_id_data= pd.read_csv("WB_id.csv")
    na_list=[]
    g_link_list=[]
    for gene in gene_multiple:
        try:
            gene_id=gene_id_data.loc[gene_id_data['Gene']==gene,'Identifier'].values[0]
            glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
            g_link_list.append(f'<a href="{glink}">{gene}</a>')
        except:
            na_list.append(gene)
    st.markdown(f"<p style='font-size:20px'>For more gene information on {', '.join(g_link_list)}(Source: GenomeAlliance)</p>", unsafe_allow_html=True)
    if na_list:
        na_links = [f'<a href="https://www.alliancegenome.org">{gene}</a>' for gene in na_list]
        st.markdown(f"<p style='font-size:20px'>Information not available for: {', '.join(na_links)}</p>", unsafe_allow_html=True)

    #filter data for particular genes
    tap_output_gene = tap_output[tap_output['Gene'].isin(gene_multiple)]
    gene_tap_data = tap_output[tap_output['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2']+ gene_multiple)].dropna(subset=['taps'])
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)
    
    #add columns for msd, habituation plots and heatmap plots
    col9, col10, col11= st.columns([1,1,1])
    tap_tstat_allele_selected = tap_tstat_allele[tap_tstat_allele['Gene'].isin(gene_multiple)]


    # Create a heatmap
    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_allele_selected.set_index('Gene').values,
        x=tap_tstat_allele_selected.set_index('Gene').columns,
        y=tap_tstat_allele_selected.set_index('Gene').index,
        colorscale='RdBu',
        zmin=-3,
        zmax=3,
        colorbar=dict(
            len=0.95,
            thickness=10,
            tickvals=[-3, 0, 3],
            ticktext=['-3', '0', '3'],
            title="",
            titleside="right"
        )
    ))

    fig.update_layout(
        width=900,
        height=1200,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title=""
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col9.subheader(f'Comprehensive heatmap of the dataset with selected genes')
    col9.plotly_chart(fig, use_container_width=True)

    # Add download buttons
    col9.download_button(
        label="Download CSV",
        data=convert_df(tap_tstat_allele_selected.set_index('Gene')),
        file_name="Data_Glance_Heatmap.csv",
        mime="text/csv",
        key='dnldheatmapcsvcustom'
    )

    col10.subheader('Rank in phenotype')
    multigene_phenotype_option = col10.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list),
        key='multigene_phenotype_select')
    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    sns.set_context('notebook')
    gene_colors = ["dimgray"] * len(gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])["Gene"])
    gene_colors[gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"] == "N2"].index[
        0]] = "red"
    for g in gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"]).reset_index(drop=True)[
        gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"]).reset_index(drop=True)["Gene"].isin(gene_multiple)].index:
        gene_colors[g] = "magenta"
    fig, ax = plt.subplots(figsize=(4, 16))
    ax = plt.errorbar(x=gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-mean"],
                    y=gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])["Gene"],
                    xerr=gene_MSD[f"{multigene_phenotype_option}-ci95_hi"] - gene_MSD[f"{multigene_phenotype_option}-mean"],
                    fmt="none", marker="none", ecolor=gene_colors, elinewidth=3)
    ax = plt.scatter(x=gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-mean"],
                    y=gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])["Gene"],
                    marker='o', color=gene_colors)
    plt.yticks(fontsize=7) # added to see the axis labels better

    plt.xlabel('Sample Mean Distance')
    plt.ylabel('Genes')
    plt.title(f"{multigene_phenotype_option}")

    multigene_phenotype_plot = io.BytesIO()
    plt.savefig(multigene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col10.image(multigene_phenotype_plot, width=None,caption=f'Sample mean distance from wildtype for selected phenotype: {multigene_phenotype_option} and selected genes :{gene_multiple}. Error bars are 95% CI.')
    
    #combine data and rename columns :
    multigene_dat=pd.concat( [gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])["Gene"],
                   gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-mean"],
                   gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-ci95_lo"],
                   gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-ci95_hi"]],
                   axis=1)
    multigene_dat.columns=["Gene", f"{multigene_phenotype_option}", f"{multigene_phenotype_option}-lower" ,f"{multigene_phenotype_option}-upper"]
    
    # Insert download graph button
    col10_1,col10_2=col10.columns(2)
    col10_1.download_button(label="Download Plot",
                        data=multigene_phenotype_plot,
                        file_name=f"multi_gene_{multigene_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldmultigenephenotypeprofile')
    col10_2.download_button(label="Download csv",
                            data=convert_df(multigene_dat),
                            file_name=f"Gene-specific Data Sample mean distance {multigene_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldmultigenephenotypeprofilecsv')
    
    col11.subheader('Habituation Curves of Response')
    genes = gene_tap_data_plot['Gene'].unique()

    # Create a list of unique colors for the genes
    colors = []
    while len(colors) < len(genes):
        color = sns.color_palette()[random.randint(1, len(sns.color_palette())-1)]
        if color not in colors:
            colors.append(color)

    # Create a palette with 'teelblue' for 'N2' and the unique colors for the other genes
    new_palette = ["steelblue" if gene == "N2" else color for gene, color in zip(genes, colors)]

    with col11:
        tab7, tab8, tab9 = st.tabs(["Probability",
                                "Duration",
                                "Speed"])
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
            tab7_1,tab7_2=tab7.columns(2)
            tab7_1.download_button(label="Download Plot",
                            data=img7,
                            file_name=f"Probability of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn13')
            tab7_2.download_button(label="Download csv",
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
            tab8_1,tab8_2=tab8.columns(2)
            tab8_1.download_button(label="Download Plot",
                            data=img8,
                            file_name=f"Duration of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn15')
            tab8_2.download_button(label="Download csv",
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
            tab9_1,tab9_2=tab9.columns(2)
            tab9_1.download_button(label="Download Plot",
                            data=img9,
                            file_name=f"Speed of Tap Habituation {gene_multiple}.png",
                            mime="image/png",
                            key='dnldbtn17')
            tab9_2.download_button(label="Download csv",
                            data=convert_df(gene_tap_data_plot),
                            file_name=f"Gene-specific Data {gene_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn18')
        # seaborn graph of Speed Habituation Curve
        # Insert download graph button
    st.download_button(label="Download raw baseline data",
                    data=convert_df(baseline_output),
                    file_name=f"raw_baseline_data.csv",
                    mime="text/csv",
                    key='dnldgenemultibaseoutcsv')


if page ==pages[4]:
   # multiple selection option for alleles
    st.header('Custom Allele Selection ')
    al_options=[]
    for a in (tap_output['dataset'].unique()):
        if a not in ['N2','cpr_5_ok2344','lfe-2-sy326']: ### for debugging this is error in backend
            al_options.append(a)
            
    allele_multiple = st.multiselect(
        label="Select Allele",
        options=al_options,
        default=al_options[0],
        placeholder="make a selection",
        help="select and de-select alleles you want to analyze",
        key="alleleselection")
    
    gene_id_data= pd.read_csv("WB_id.csv")# data for gene id
    allele_id_data= pd.read_csv("Gene_Allele_WormBaseID.csv",names=['Identifier', 'Gene', 'Allele'])# data for allele id
    na_list=[]
    g_link_list=[] # list for gene links (AllianceGenome)
    w_link_list=[] # list for allele links (WormBase)
    allele_list=[]
    for a in allele_multiple:
        gene, allele= a.split('_')
        allele_list.append(a)
        gene_id = gene_id_data.loc[gene_id_data['Gene'] == gene, 'Identifier']
        allele_id = allele_id_data.loc[(allele_id_data['Gene'] == gene) & (allele_id_data['Allele'] == allele), 'Identifier']
        
        if allele_id.any():
            allele_id=allele_id.values[0]
            wlink=f'https://wormbase.org/species/c_elegans/variation/{allele_id}'
            w_link_list.append(f'<a href="{wlink}">{allele}</a>')

        else: # if allele option doesnt match then send to default page
            allele_id= "default value"
            na_list.append(allele)
        
        if gene_id.any():
            gene_id=gene_id.values[0]
            glink=f'https://www.alliancegenome.org/gene/WB:{gene_id}'
            g_link_list.append(f'<a href="{glink}">{gene}</a>')

        else:
            gene_id= "default value"
            na_list.append(gene)

    st.markdown(f"<p style='font-size:20px'>For more gene information on {', '.join(g_link_list)}(Source: GenomeAlliance) </p>", unsafe_allow_html=True)
    st.markdown(f"<p style='font-size:20px'>For more allele information on {', '.join(w_link_list)}(Source: WormBase) </p>", unsafe_allow_html=True)
    if na_list:
        na_links = [f'<a href="https://www.alliancegenome.org">{gene}</a>' for gene in na_list]
        na_links = [f'<a href="https://wormbase.org/species/c_elegans/variation/">{allele}</a>' for allele in na_list]

        st.markdown(f"<p style='font-size:20px'>Information not available for: {', '.join(na_links)}</p>", unsafe_allow_html=True)

    #filter data for particular allele

    tap_output_allele = tap_output[tap_output['dataset'].isin(allele_multiple)]
    allele_tap_data = tap_output[tap_output['Date'].isin(tap_output_allele['Date'].unique())]
    allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2']+ allele_multiple)].dropna(subset=['taps'])
    allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)
    # st.write(allele_tap_data_plot)
  
    #add columns for msd, habituation plots and heatmap plots
    col12, col13, col14= st.columns([1,1,1])
    # Filter the dataframe for the selected genes
    tap_tstat_data_selected = tap_tstat_data[tap_tstat_data['dataset'].isin(allele_list)]
    
    # Create a heatmap
    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_data_selected.set_index('dataset').values,
        x=tap_tstat_data_selected.set_index('dataset').columns,
        y=tap_tstat_data_selected.set_index('dataset').index,
        colorscale='RdBu',
        zmin=-3,
        zmax=3,
        colorbar=dict(
            len=0.95,
            thickness=10,
            tickvals=[-3, 0, 3],
            ticktext=['-3', '0', '3'],
            title="",
            titleside="right"
        )
    ))

    fig.update_layout(
        width=900,
        height=1200,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title=""
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col12.subheader(f'Comprehensive heatmap of the dataset with selected alleles')
    col12.plotly_chart(fig, use_container_width=True)

    # Add download buttons
    col12.download_button(
        label="Download CSV",
        data=convert_df(tap_tstat_data_selected.set_index('dataset')),
        file_name="Data_Glance_Heatmap.csv",
        mime="text/csv",
        key='dnldheatmapcsvcustomallele'
    )

    col13.subheader('Rank in phenotype')
    multiallele_phenotype_option = col13.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list),
        key='multiallele_phenotype_select')
    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    sns.set_context('notebook')
    allele_colors = ["dimgray"] * len(allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])["dataset"])
    allele_colors[allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"]).reset_index(drop=True)[
        allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"]).reset_index(drop=True)["dataset"] == "N2"].index[
        0]] = "red"
    for a in allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"]).reset_index(drop=True)[
        allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"]).reset_index(drop=True)["dataset"].isin(allele_multiple)].index:
        allele_colors[a] = "magenta"
    fig, ax = plt.subplots(figsize=(4, 16))
   
    ax = plt.errorbar(x=allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-mean"],
                    y=allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])["dataset"],
                    xerr=allele_MSD[f"{multiallele_phenotype_option}-ci95_hi"] - allele_MSD[f"{multiallele_phenotype_option}-mean"],
                    fmt="none", marker="none", ecolor=allele_colors, elinewidth=3)
    ax = plt.scatter(x=allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-mean"],
                    y=allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])["dataset"],
                    marker='o', color=allele_colors)
    plt.yticks(fontsize=7) # added to see the axis labels better

    plt.xlabel('Sample Mean Distance')
    plt.ylabel('Gene_Allele')
    plt.title(f"{multiallele_phenotype_option}")
    
    #edit from here
    multiallele_phenotype_plot = io.BytesIO()
    plt.savefig(multiallele_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    #display image 
    col13.image(multiallele_phenotype_plot, width=None,caption=f'Sample mean distance from wildtype for selected phenotype: {multiallele_phenotype_option} and selected alleles :{allele_multiple}. Error bars are 95% CI.')
    
    #combine data and rename columns :
    multiallele_dat=pd.concat( [allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])["dataset"],
                   allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-mean"],
                   allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-ci95_lo"],
                   allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-ci95_hi"]],
                   axis=1)
    multiallele_dat.columns=["Allele", f"{multiallele_phenotype_option}", f"{multiallele_phenotype_option}-lower" ,f"{multiallele_phenotype_option}-upper"]
    
    # Insert download graph button
    col13_1,col13_2 = col13.columns(2)
    col13_1.download_button(label="Download Plot",
                        data=multiallele_phenotype_plot,
                        file_name=f"multi_allele_{multiallele_phenotype_option}_profile.png",
                        mime="image/png",
                        key='dnldmultiallelephenotypeprofile')
    col13_2.download_button(label="Download csv",
                            data=convert_df(multiallele_dat),
                            file_name=f"Allele-specific Data Sample mean distance {multiallele_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldmultiallelephenotypeprofilecsv')
    
    col14.subheader('Habituation Curves of Response')
    alleles = allele_tap_data_plot['dataset'].unique()

    # Create a cycle of unique colors
    colors_list = sns.color_palette()[1:] # excludes steelblue from palette
    color_cycle = itertools.cycle(colors_list)

    # Create a list of unique colors for the alleles
    colors = [next(color_cycle) for _ in range(len(alleles))]

    # Create a palette with 'teelblue' for 'N2' and the unique colors for the other genes
    new_palette = ["steelblue" if allele == "N2" else color for allele, color in zip(alleles, colors)]

    with col14:
        tab10, tab11, tab12 = st.tabs(["Probability",
                                "Duration",
                                "Speed"])
        with tab10:
            #  Habituation of Response Probability Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            plt.gca().xaxis.grid(False)  # <- gets rid of x-axis markers to make data look clean
            ax = sns.pointplot(x="taps",  # <- Here we use seaborn as our graphing package.
                            y="prob",
                            data=allele_tap_data_plot,
                            hue='Strain',  # <- Here we use the extra column from step 6 to separate by group
                            palette=new_palette,
                            errorbar='se')  # <- Confidence interval. 95 = standard error
            plt.xlabel("Taps")  # <- x-axis title
            plt.ylabel("Probability")  # <- y-axis title
            plt.title("Habituation of Response Probability", fontsize='16')  # <- Figure Title
            plt.ylim(0, 1)
            ax.legend(loc='upper right', fontsize='12')  # <- location of your legend

            # download graph button
            img10 = io.BytesIO()
            plt.savefig(img10, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab10.image(img10, width=None,caption=(f'Habituation of Response Probability: {allele_multiple}'))
            # Insert download plot and download csv button
            tab10_1,tab10_2=tab10.columns(2)
            tab10_1.download_button(label="Download Plot",
                            data=img10,
                            file_name=f"Probability of Tap Habituation {allele_multiple}.png",
                            mime="image/png",
                            key='dnldbtn19')
            tab10_2.download_button(label="Download csv",
                            data=convert_df(allele_tap_data_plot),
                            file_name=f"Allele-specific Data {allele_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn20')

        with tab11:
            #  Habituation of Response Duration Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="dura",
                            data=allele_tap_data_plot,
                            hue='Strain',
                            palette=new_palette, # N2 to be blue consistently
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Duration", fontsize='12')
            plt.title("Habituation of Response Duration", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
            
            # download graph button
            img11 = io.BytesIO()
            plt.savefig(img11, format='png', dpi=300, bbox_inches='tight')
            #display image 
            tab11.image(img11, width=None,caption=(f'Habituation of Response Duration: {allele_multiple}'))
            # Insert download plot and download csv button
            tab11_1,tab11_2=tab11.columns(2)
            tab11_1.download_button(label="Download Plot",
                            data=img11,
                            file_name=f"Duration of Tap Habituation {allele_multiple}.png",
                            mime="image/png",
                            key='dnldbtn21')
            tab11_2.download_button(label="Download csv",
                            data=convert_df(allele_tap_data_plot),
                            file_name=f"Allele-specific Data {allele_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn22')

        with tab12:
            #  Habituation of Response Speed Plot
            fig, ax = plt.subplots(figsize=(12, 10))
            # seaborn plot
            ax = sns.pointplot(x="taps",
                            y="speed",
                            data=allele_tap_data_plot,
                            hue='Strain',
                            palette=new_palette, # N2 to be blue consistently
                            errorbar='se')
            plt.xlabel("Taps", fontsize='12')
            plt.ylabel("Speed", fontsize='12')
            plt.title("Habituation of Response Speed", fontsize='16')
            plt.ylim(0, None)
            ax.legend(loc='upper right', fontsize='12')
        
            img12 = io.BytesIO()
            plt.savefig(img12, format='png', dpi=300, bbox_inches='tight')
            
            #display image 
            tab12.image(img12, width=None,caption=(f'Habituation of Response Speed: {allele_multiple}'))
            
            # Insert download plot and download csv button
            tab12_1,tab12_2=tab12.columns(2)
            tab12_1.download_button(label="Download Plot",
                            data=img12,
                            file_name=f"Speed of Tap Habituation {allele_multiple}.png",
                            mime="image/png",
                            key='dnldbtn23')
            tab12_2.download_button(label="Download csv",
                            data=convert_df(allele_tap_data_plot),
                            file_name=f"Allele-specific Data {allele_multiple}.csv",
                            mime="text/csv",
                            key='dnldbtn24')
    st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldallelemultibaseoutcsv')
    
        

    