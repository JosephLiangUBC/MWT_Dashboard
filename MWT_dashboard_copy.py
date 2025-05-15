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
import hmac
import psycopg

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if hmac.compare_digest(st.session_state["password"], st.secrets["password"]):
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # Don't store the password.
        else:
            st.session_state["password_correct"] = False

    # Return True if the password is validated.
    if st.session_state.get("password_correct", False):
        return True

    # Show input for password.
    st.text_input(
        "Password", type="password", on_change=password_entered, key="password"
    )
    if "password_correct" in st.session_state:
        st.error("😕 Password incorrect")
    return False


if not check_password():
    st.stop()  # Do not continue if check_password is not True.


# Added configuration to wide for desktop screen (do not remove)
# it has to be the first command in the file for it to run
# Streamlit command starts here
st.set_page_config(layout="wide")

# convert dataframe to csv for download
@st.cache_data
def convert_df(df):
    return df.to_csv().encode("utf-8")

# Fail gracefully option
def read(table, connection):
    with connection.cursor() as cursor:
        cursor.execute(f'SELECT * from "{table}"')

        # Fetch all rows from database
        record = cursor.fetchall()
        column_names = [desc[0] for desc in cursor.description]
    return pd.DataFrame(data=record, columns=column_names)

@st.cache_data
def aggregate_unique_values(df ,by):
    """Aggregate and transform tstat_gene_data, tstat_allele_data , gene_profile_data, and  allele_profile_data table"""
    if(len(by)>1):    
        df['Metric'] = pd.Categorical(df['Metric'], categories=df['Metric'].unique(), ordered=True)

    grouped = df.groupby(by).agg(
        {col: 'mean' for col in df.columns if col not in by+ ['Screen']}).reset_index()

    # Aggregate the Screen column into a list
    grouped['Screen'] = df.groupby(by)['Screen'].apply(lambda x: list(set(x))).reset_index(drop=True)

    return grouped

@st.cache_data
def aggregate_unique_values_MSD(df, by):
    """Aggregate and transform gene_MSD and allele_MSD table"""
    # Define the columns to aggregate
    agg_cols = [col for col in df.columns if col not in by + ['Screen']]

    # Aggregate the columns using a weighted average
    grouped = df.groupby(by).apply(lambda x: pd.Series({
        col: (x[col] * x[col.replace('-mean', '-count')]).sum() / x[col.replace('-mean', '-count')].sum()
        if '-mean' in col
        else np.sqrt((x[col] ** 2 * x[col.replace('-sem', '-count')]).sum() / x[col.replace('-sem', '-count')].sum())
        if '-sem' in col
        else x[col].sum()
        if '-count' in col
        else np.nan
        for col in agg_cols
    })).reset_index()

    # Calculate new confidence intervals
    for col in grouped.columns:
        if '-mean' in col:
            mean_col = col
            sem_col = col.replace('-mean', '-sem')
            count_col = col.replace('-mean', '-count')
            ci95_lo_col = col.replace('-mean', '-ci95_lo')
            ci95_hi_col = col.replace('-mean', '-ci95_hi')

            grouped[ci95_lo_col] = grouped[mean_col] - 1.96 * grouped[sem_col] / np.sqrt(grouped[count_col])
            grouped[ci95_hi_col] = grouped[mean_col] + 1.96 * grouped[sem_col] / np.sqrt(grouped[count_col])

    # Aggregate the Screen column into a list
    grouped['Screen'] = df.groupby(by)['Screen'].apply(lambda x: list(set(x))).reset_index(drop=True)

    return grouped

@st.cache_data
def fetch_data():
    with psycopg.connect(dbname="mwtdata", user=st.secrets["psql_user"], password=st.secrets["psql_passwword"], host="rds-mwt-data.ctie02ksmcqc.ca-central-1.rds.amazonaws.com", port=5432) as connection:
        # Read data from SQLite database
        tap_output = read('tap_response_data', connection)
        tap_tstat_allele = aggregate_unique_values(read('tstat_gene_data', connection),["Gene"]).explode('Screen').reset_index(drop=True)
        tap_tstat_data = aggregate_unique_values(read('tstat_allele_data', connection),["dataset"]).explode('Screen').reset_index(drop=True)
        # allele_metric_data = read('allele_phenotype_data')
        gene_profile_data = aggregate_unique_values(read('gene_profile_data', connection),['Gene','Metric']).explode('Screen').reset_index(drop=True)
        allele_profile_data = aggregate_unique_values(read('allele_profile_data', connection),['dataset','Metric']).explode('Screen').reset_index(drop=True)
        gene_MSD = aggregate_unique_values_MSD(read('gene_MSD', connection),["Gene"]).explode('Screen').reset_index(drop=True)
        allele_MSD = aggregate_unique_values_MSD(read('allele_MSD', connection),["dataset"]).explode('Screen').reset_index(drop=True)
        id_data=read('Gene_Allele_WormBaseID', connection) ##table in database with wormbase id's for all genes and alleles
        
    return tap_output, tap_tstat_allele, tap_tstat_data, gene_profile_data, allele_profile_data, gene_MSD, allele_MSD, id_data

tap_output, tap_tstat_allele, tap_tstat_data, gene_profile_data, allele_profile_data, gene_MSD, allele_MSD, id_data = fetch_data()
tap_output['Strain'] = tap_output['Gene'] + " (" + tap_output['Allele'] + ")"

# Defining the color palette for the plots
metric_palette = ["k", "k", "k",
                  "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgrey","darkgray",
                  "lightsteelblue", "lightsteelblue", "lightsteelblue",
                  "powderblue", "powderblue", "powderblue",
                  "cadetblue", "cadetblue", "cadetblue",
                  "thistle", "thistle", "thistle",
                  "slateblue","slateblue","slateblue"]

# config setting for plotly
config = {
  'toImageButtonOptions': {
    'format': 'png', # one of png, svg, jpeg, webp
    'filename': 'Image',
    'height': None,
    'width': None,
    'scale': 3 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

# creating tabs for dashboard
pages = ["Home - Data at a Glance", "Gene-specific Data", "Allele-specific Data",  "Custom Gene Selection","Custom Allele Selection", "Help & Documentation", "Citations"]
# to-do: add 'clustering' in the pages list above at 5th position 
page = st.sidebar.radio("Select a page", pages)

# Streamlit Dashboard title
st.title('MWT Data Dashboard - Rankin Lab @ UBC')


def select_datasets():
    global gene_MSD
    global phenotype_list
    global tap_output
    global tap_tstat_allele
    global tap_tstat_data
    global gene_profile_data
    global allele_profile_data
    global allele_MSD

    # Select dataset option
    datasets = st.multiselect(
        label="Select Datasets",
        options=gene_MSD.Screen.unique(),
        default="PD_Screen",
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
    dropna_features.remove('Memory Retention of Response Duration')
    dropna_features.remove('Memory Retention of Response Probability')
    dropna_features.remove('Memory Retention of Response Speed')
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

# Visualisations for data tab
if page ==pages[0]:
    st.header('Home - Data at a Glance')
    select_datasets()

    col1, col2 = st.columns([4, 5])

    col1.subheader("For A Single Phenotype")

    phenotype_option = col1.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list), key="phenotypeselect")



    data_sorted = gene_MSD.sort_values(by=[f"{phenotype_option}-mean"]).reset_index(drop=True)
    colors = ["dimgrey"] * len(data_sorted)

    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        # color = colors[i]
        if row['Gene'] == "N2":
            colors="red"
        else:
            colors="dimgray"
        fig.add_trace(go.Scatter(
            x=[row[f"{phenotype_option}-mean"]],
            y=[row["Gene"]],
            error_x=dict(
                type='data',
                array=[row[f"{phenotype_option}-ci95_hi"] - row[f"{phenotype_option}-mean"]],
                arrayminus=[row[f"{phenotype_option}-mean"] - row[f"{phenotype_option}-ci95_lo"]],
                visible=True,
                color=colors,
                thickness=3,
                width=0
            ),
            mode='markers',
            marker=dict(
                color=colors,
                size=12,
                symbol='circle',
                line=dict(
                    color='rgb(0,0,0)',
                    width=1
                ),
            ),
            showlegend=False,  # Hide individual points from legend
            name=""
            
        ))
    # Add vertical line at 0
    fig.add_vline(x=0,  line_width=1, line_dash="dash", line_color="red")
    # Update layout with labels and title
    fig.update_layout(
        title={
            'text': f"{phenotype_option}"},
        xaxis_title='Sample Mean Distance',
        yaxis_title='Gene',
        plot_bgcolor='white',
        paper_bgcolor='white',
        font_color="black",
        width=600,
        height=800,
        yaxis=dict(showticklabels=True, 
                   dtick=1,
                   tickfont=dict(color='black', size=6)),
        margin=dict(l=10, r=10, t=40, b=0),  # Adjust margins as needed
        annotations=[
            dict(
                text=f'Sample mean distance from wildtype for all strains for selected phenotype: {phenotype_option}. Error bars are 95% CI',
                xref="paper",
                yref="paper",
                x=0,
                y=-0.2,
                showarrow=False,
                font=dict(
                    size=12,
                    color="black"
                )
            )
        ]
    )
    
    phenotype_plot = io.BytesIO()
    fig.write_image(phenotype_plot, format='png',scale=3)
    phenotype_plot.seek(0)
    col1.plotly_chart(fig, use_container_width=True, **{'config': config})
    
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
                            file_name=f"Data Glance Sample Mean Distance {phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldphenotypeprofilecsv')

    # Insert download graph button

    # Create a heatmap
    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_allele.set_index('Gene').values,
        x=tap_tstat_allele.set_index('Gene').columns,
        y=tap_tstat_allele.set_index('Gene').index,
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
        yaxis_title="",
        yaxis=dict(showticklabels=True, 
                   dtick=1,
                   tickfont=dict(color='black', size=6)),
        xaxis=dict(showticklabels=True, 
                   dtick=1,
                   tickfont=dict(color='black', size=12))
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col2.subheader("Comprehensive Heatmap of Entire Dataset")
    col2.plotly_chart(fig, use_container_width=True, **{'config': config})

    col2_1,col2_2= col2.columns(2)
    col2_1.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmap')
    # Add download buttons    
    col2_2.download_button(
        label="Download CSV",
        data=convert_df(tap_tstat_allele.set_index('Gene')),
        file_name="Data_Glance_Heatmap.csv",
        mime="text/csv",
        key='dnldheatmapcsv'
    )

    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readbaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
        
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldbaseoutcsv')

if page == pages[1]:
    st.header('Gene-specific Data')
    select_datasets()
    # Create a session state for the gene selection
    st.session_state.setdefault('gene_select', None)
    gene_option = st.selectbox(
        'Select a gene',
        [gene for gene in tap_output['Gene'].unique()if gene != 'N2'], 
        key="geneselect",
        index=[gene for gene in tap_output['Gene'].unique() if gene != 'N2'].index(st.session_state.gene_select) if st.session_state.gene_select else 0        
        )
    
    st.session_state.gene_select = gene_option

    if gene_option:
        gene_id = id_data.loc[id_data['Gene'] == gene_option, 'WBGene'].values
        if len(gene_id) == 0:
            gene_id = id_data.loc[id_data['Sequence'] == gene_option, 'WBGene'].values

        if len(gene_id) > 0:
            glink = f'https://www.alliancegenome.org/gene/WB:{gene_id[0]}'
            st.markdown(f'<p style="font-size:20px">For more gene information on <a href="{glink}">{gene_option}</a> (Source: GenomeAlliance)</p>', unsafe_allow_html=True)
        else:
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
    data_sorted = gene_MSD.sort_values(by=[f"{gene_phenotype_option}-mean"])
    gene_colors = ["dimgray"] * len(data_sorted["Gene"])
    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            gene_colors="red"
        elif row['Gene'] == gene_option:
            gene_colors="magenta"
        else:
            gene_colors="dimgray"
        # color = "red" if row['Gene'] == "N2" else "dimgrey"
        fig.add_trace(go.Scatter(
            x=[row[f"{gene_phenotype_option}-mean"]],
            y=[row["Gene"]],
            error_x=dict(
                type='data',
                array=[row[f"{gene_phenotype_option}-ci95_hi"] - row[f"{gene_phenotype_option}-mean"]],
                arrayminus=[row[f"{gene_phenotype_option}-mean"] - row[f"{gene_phenotype_option}-ci95_lo"]],
                visible=True,
                color=gene_colors,
                thickness=3,
                width=0
            ),
            mode='markers',
            marker=dict(
                color=gene_colors,
                size=12,
                symbol='circle',
                line=dict(
                    color='rgb(0,0,0)',
                    width=1
                ),
            ),
            showlegend=False,  # Hide individual points from legend
            name=""
            
        ))
    # Add vertical line at 0
    fig.add_vline(x=0,  line_width=1, line_dash="dash", line_color="red")
    # Update layout with labels and title
    fig.update_layout(
        title=f"{gene_phenotype_option}",
        xaxis_title='Sample Mean Distance',
        yaxis_title='Gene',
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=600,
        height=1200,
        yaxis=dict(showticklabels=True, 
                   dtick=1,
                   tickfont=dict(color='black', size=6)),
        margin=dict(l=100, r=50, t=100, b=50),  # Adjust margins as needed
        annotations=[
            dict(
                text=f'Sample mean distance from wildtype for all strains for selected phenotype: {gene_phenotype_option}. Error bars are 95% CI',
                xref="paper",
                yref="paper",
                x=0,
                y=-0.2,
                showarrow=False,
                font=dict(
                    size=12,
                    color="black"
                )
            )
        ]
    )


    gene_phenotype_plot = io.BytesIO()
    fig.write_image(gene_phenotype_plot, format='png',scale=3)
    gene_phenotype_plot.seek(0)
    col4.plotly_chart(fig, use_container_width=True, **{'config': config})

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
                            palette=["black" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
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
                            palette=["black" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
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
                            palette=["black" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
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

    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readgenebaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        for conn_path in conn_list:
            try:
                conn = sqlite3.connect(conn_path) 
                break
            except:
                pass
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldgenebaseoutcsv')

if page ==pages[2]:
    st.header('Allele-specific Data')
    select_datasets()
    # select allele 
    st.session_state.setdefault('allele_select', None)
    allele_option = st.selectbox(
        'Select a allele',
        [allele for allele in tap_output['dataset'].unique() if allele != 'N2'],
        index=[allele for allele in tap_output['dataset'].unique() if allele != 'N2'].index(st.session_state.allele_select) if st.session_state.allele_select else 0
        )
    
    st.session_state.allele_select = allele_option

    #splititing gene allele to get gene and allele columns
    gene, allele = allele_option.split('_')

    # check if allele and gene options match the preexisting list
    if allele_option:
        gene_id = id_data.loc[id_data['Gene'] == gene, 'WBGene']
        allele_id = id_data.loc[(id_data['Gene'] == gene) & (id_data['Allele'] == allele), 'WBAllele']
        allele_id = (allele_id.values[0]) if allele_id.any() else "default value"
        gene_id=(gene_id.values[0]) if gene_id.any() else "default value"


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
    
    data_sorted = allele_MSD.sort_values(by=[f"{allele_phenotype_option}-mean"])
    allele_colors = ["dimgray"] * len(data_sorted["dataset"])
    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['dataset'] == "N2":
            allele_colors="red"
        elif row['dataset'] == allele_option:
            allele_colors="magenta"
        else:
            allele_colors="dimgrey"
        # color = "red" if row['Gene'] == "N2" else "dimgrey"
        fig.add_trace(go.Scatter(
            x=[row[f"{allele_phenotype_option}-mean"]],
            y=[row["dataset"]],
            error_x=dict(
                type='data',
                array=[row[f"{allele_phenotype_option}-ci95_hi"] - row[f"{allele_phenotype_option}-mean"]],
                arrayminus=[row[f"{allele_phenotype_option}-mean"] - row[f"{allele_phenotype_option}-ci95_lo"]],
                visible=True,
                color=allele_colors,
                thickness=3,
                width=0
            ),
            mode='markers',
            marker=dict(
                color=allele_colors,
                size=12,
                symbol='circle',
                line=dict(
                    color='rgb(0,0,0)',
                    width=1
                ),
            ),
            name="",
            showlegend=False  # Hide individual points from legend
        ))
    # Add vertical line at 0
    fig.add_vline(x=0,  line_width=1, line_dash="dash", line_color="red")
    # Update layout with labels and title
    fig.update_layout(
        title=f"{allele_phenotype_option}",
        xaxis_title='Sample Mean Distance',
        yaxis_title='Gene',
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=600,
        height=1200,
        yaxis=dict(showticklabels=True, 
                   dtick=1,
                   tickfont=dict(color='black', size=6)),
        margin=dict(l=100, r=50, t=100, b=50),  # Adjust margins as needed
        annotations=[
            dict(
                text=f'Sample mean distance from wildtype for all strains for selected phenotype: {allele_phenotype_option}. Error bars are 95% CI',
                xref="paper",
                yref="paper",
                x=0,
                y=-0.2,
                showarrow=False,
                font=dict(
                    size=12,
                    color="black"
                )
            )
        ]
    )

    allele_phenotype_plot = io.BytesIO()
    fig.write_image(allele_phenotype_plot, format='png',scale=3)
    allele_phenotype_plot.seek(0)
    col6.plotly_chart(fig, use_container_width=True, **{'config': config})

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
                            palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
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
                            palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
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
                            palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
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


    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readallelebaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        for conn_path in conn_list:
            try:
                conn = sqlite3.connect(conn_path) 
                break
            except:
                pass
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldallelebaseoutcsv')

if page ==pages[3]:
   # multiple selection option for genes
    st.header('Custom Gene Selection ')
    select_datasets()
    st.session_state.setdefault('gene_select', [gene for gene in tap_output['Gene'].unique() if gene != 'N2'][0])

    gene_multiple = st.multiselect(
        label="Select Genes",
        options=[gene for gene in tap_output['Gene'].unique() if gene != 'N2'],
        default=st.session_state.gene_select,
        placeholder="make a selection",
        help="select and de-select genes you want to analyze",
        key="geneselection")
    st.session_state.gene_select = gene_multiple

    na_list=[]
    g_link_list=[]
    for gene in gene_multiple:
        gene_id = id_data.loc[id_data['Gene'] == gene, 'WBGene'].values
        if len(gene_id) == 0:
            gene_id = id_data.loc[id_data['Sequence'] == gene, 'WBGene'].values
        if len(gene_id) > 0:
            glink = f'https://www.alliancegenome.org/gene/WB:{gene_id[0]}'
            g_link_list.append(f'<a href="{glink}">{gene}</a>')
        else:
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
    h = 200*len(gene_multiple)
    fig.update_layout(
        width=900,
        height=h,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(tickangle=0),
        xaxis=dict(showticklabels=True,tickfont=dict(size=8))
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col9.subheader(f'Comprehensive heatmap of the dataset with selected genes')
    col9.plotly_chart(fig, use_container_width=True, **{'config': config})

    # Add download buttons
    col9_1,col9_2=col9.columns(2)
    col9_1.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmapcustom')
    col9_2.download_button(
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
    data_sorted=gene_MSD.sort_values(by=[f"{multigene_phenotype_option}-mean"])
    gene_colors = ["dimgray"] * len(data_sorted["Gene"])
    ticktext = []
    tickvals = []
    fig = go.Figure()
    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            gene_colors="red"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        elif row['Gene'] in (gene_multiple):
            gene_colors="magenta"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        else:
            gene_colors="dimgray"
        # color = "red" if row['Gene'] == "N2" else "dimgrey"
        fig.add_trace(go.Scatter(
            x=[row[f"{multigene_phenotype_option}-mean"]],
            y=[row["Gene"]],
            error_x=dict(
                type='data',
                array=[row[f"{multigene_phenotype_option}-ci95_hi"] - row[f"{multigene_phenotype_option}-mean"]],
                arrayminus=[row[f"{multigene_phenotype_option}-mean"] - row[f"{multigene_phenotype_option}-ci95_lo"]],
                visible=True,
                color=gene_colors,
                thickness=3,
                width=0
            ),
            mode='markers',
            marker=dict(
                color=gene_colors,
                size=12,
                symbol='circle',
                line=dict(
                    color='rgb(0,0,0)',
                    width=1
                ),
            ),
            showlegend=False,  # Hide individual points from legend
            name=""
            
        ))

        # Update layout with labels and title
        fig.update_layout(
            title=f"{multigene_phenotype_option}",
            xaxis_title='Sample Mean Distance',
            yaxis_title='Gene',
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=600,
            height=1200,
            yaxis=dict(
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext),
            margin=dict(l=100, r=50, t=100, b=50),  # Adjust margins as needed
            annotations=[
                dict(
                    text=f'Sample mean distance from wildtype for all strains for selected phenotypes: {multigene_phenotype_option}. Error bars are 95% CI',
                    xref="paper",
                    yref="paper",
                    x=0,
                    y=-0.2,
                    showarrow=False,
                    font=dict(
                        size=12,
                        color="black"
                    )
                )
            ]
        )
    multigene_phenotype_plot = io.BytesIO()
    fig.write_image(multigene_phenotype_plot, format='png',scale=3)
    multigene_phenotype_plot.seek(0)
    col10.plotly_chart(fig, use_container_width=True, **{'config': config})
   
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
                            data=convert_df(multigene_dat[multigene_dat['Gene'].isin(gene_multiple)]),
                            file_name=f"Gene-specific Data Sample mean distance {multigene_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldmultigenephenotypeprofilecsv')
    
    col11.subheader('Habituation Curves of Response')
    genes = gene_tap_data_plot['Gene'].unique()
    # Create a cycle of unique colors
    colors_list = sns.color_palette("husl",n_colors=len(genes)+1)[:] 
    color_cycle = itertools.cycle(colors_list)

    # Create a list of unique colors for the alleles
    colors = [next(color_cycle) for _ in range(len(genes))]

    # Create a palette with 'black' for 'N2' and the unique colors for the other genes
    new_palette = ["black" if gene == "N2" else color for gene, color in zip(genes, colors)]

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


    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readgenemultibaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        for conn_path in conn_list:
            try:
                conn = sqlite3.connect(conn_path) 
                break
            except:
                pass
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output[baseline_output['Gene'].isin(gene_multiple)]),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldgenemultibaseoutcsv')


if page ==pages[4]:
   # multiple selection option for alleles
    st.header('Custom Allele Selection')
    select_datasets()
    st.session_state.setdefault('allele_select', [allele for allele in tap_output['dataset'].unique() if allele != 'N2'][0])

    allele_multiple = st.multiselect(
        label="Select Allele",
        options=[allele for allele in tap_output['dataset'].unique() if allele != 'N2'],
        default=st.session_state.allele_select,
        placeholder="make a selection",
        help="select and de-select alleles you want to analyze",
        key="alleleselection")
    st.session_state.allele_select = allele_multiple

    na_list=[]
    g_link_list=[] # list for gene links (AllianceGenome)
    w_link_list=[] # list for allele links (WormBase)
    allele_list=[]
    for a in allele_multiple:
        gene, allele= a.split('_')
        allele_list.append(a)
        gene_id = id_data.loc[id_data['Gene'] == gene, 'WBGene']
        allele_id = id_data.loc[(id_data['Gene'] == gene) & (id_data['Allele'] == allele), 'WBAllele']
        
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
    h = 200*len(allele_multiple)
    fig.update_layout(
        width=900,
        height=h,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(tickangle=0),
        xaxis=dict(showticklabels=True,tickfont=dict(size=8))
    )

    imgheatmap = io.BytesIO()
    #imgheatmap=fig.to_image(imgheatmap, format='png', width=1100, height=h+200, scale=3)
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    # Display the heatmap in Streamlit
    col12.subheader(f'Comprehensive heatmap of the dataset with selected alleles')
    col12.plotly_chart(fig, use_container_width=True, **{'config': config})

    # Add download buttons
    col12_1,col12_2=col12.columns(2)
    col12_1.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmapcustomallele')
    col12_2.download_button(
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
    data_sorted= allele_MSD.sort_values(by=[f"{multiallele_phenotype_option}-mean"])
    allele_colors = ["dimgray"] * len(data_sorted["dataset"])
    ticktext = []
    tickvals = []
    fig = go.Figure()
    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['dataset'] == "N2":
            allele_colors="red"
            ticktext.append(row['dataset'])
            tickvals.append(row['dataset'])
        elif row['dataset'] in (allele_multiple):
            allele_colors="magenta"
            ticktext.append(row['dataset'])
            tickvals.append(row['dataset'])
        else:
            allele_colors="dimgray"
        # color = "red" if row['Gene'] == "N2" else "dimgrey"
        fig.add_trace(go.Scatter(
            x=[row[f"{multiallele_phenotype_option}-mean"]],
            y=[row["dataset"]],
            error_x=dict(
                type='data',
                array=[row[f"{multiallele_phenotype_option}-ci95_hi"] - row[f"{multiallele_phenotype_option}-mean"]],
                arrayminus=[row[f"{multiallele_phenotype_option}-mean"] - row[f"{multiallele_phenotype_option}-ci95_lo"]],
                visible=True,
                color=allele_colors,
                thickness=3,
                width=0
            ),
            mode='markers',
            marker=dict(
                color=allele_colors,
                size=12,
                symbol='circle',
                line=dict(
                    color='rgb(0,0,0)',
                    width=1
                ),
            ),
            showlegend=False,  # Hide individual points from legend
            name=""
            
        ))

        # Update layout with labels and title
        fig.update_layout(
            title=f"{multiallele_phenotype_option}",
            xaxis_title='Sample Mean Distance',
            yaxis_title='Gene',
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=600,
            height=1200,
            yaxis=dict(
            tickmode='array',
            tickvals=tickvals,
            ticktext=ticktext),
            margin=dict(l=100, r=50, t=100, b=50),  # Adjust margins as needed
            annotations=[
                dict(
                    text=f'Sample mean distance from wildtype for all strains for selected phenotypes: {multiallele_phenotype_option}. Error bars are 95% CI',
                    xref="paper",
                    yref="paper",
                    x=0,
                    y=-0.2,
                    showarrow=False,
                    font=dict(
                        size=12,
                        color="black"
                    )
                )
            ]
        )
    
    multiallele_phenotype_plot = io.BytesIO()
    fig.write_image(multiallele_phenotype_plot, format='png',scale=3)
    multiallele_phenotype_plot.seek(0)
    col13.plotly_chart(fig, use_container_width=True, **{'config': config})
   
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
                            data=convert_df(multiallele_dat[multiallele_dat['Allele'].isin(allele_list)]),
                            file_name=f"Allele-specific Data Sample mean distance {multiallele_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldmultiallelephenotypeprofilecsv')
    
    col14.subheader('Habituation Curves of Response')
    alleles = allele_tap_data_plot['dataset'].unique()

    # Create a cycle of unique colors
    colors_list = sns.color_palette("husl",n_colors=len(alleles)+1)[:] # excludes steelblue from palette
    color_cycle = itertools.cycle(colors_list)

    # Create a list of unique colors for the alleles
    colors = [next(color_cycle) for _ in range(len(alleles))]

    # Create a palette with 'black' for 'N2' and the unique colors for the other genes
    new_palette = ["black" if allele == "N2" else color for allele, color in zip(alleles, colors)]

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
                            hue='dataset',  # <- Here we use the extra column from step 6 to separate by group
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
                            hue='dataset',
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
                            hue='dataset',
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
            

    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readallelemultibaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        for conn_path in conn_list:
            try:
                conn = sqlite3.connect(conn_path) 
                break
            except:
                pass
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(datasets)].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output[baseline_output['dataset'].isin(allele_multiple)]),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldallelemultibaseoutcsv')
    
if page ==pages[5]:
    st.markdown("""
    ## Help and Documentation

    ### About the Data Dashboard:
    Since the inception of the [Multi-Worm Tracke (MWT)](https://doi.org/10.1038/nmeth.1625), the Rankin Lab has collected vast amounts of data from a large number of strains in a number of screens that have been conducted by members of the lab. 
    This dashboard is a tool written by Joseph Liang (PhD Candidate) in an attempt to consolidate the MWT data collected over the years that follow the standard 10-second ISI (30 stimuli at 10-second inter-stimulus intervals) habituation protocols employed by members in the lab.
    Although the MWT is very powerful in its ability to gather rich multi-phenotype data from large populations of living animals simultaneously, raw data from the MWT tracker is generally not very accessible to users not familiar with the innor workings of the MWT.
    This tool seeks to solve this problem by consolidating all the appropriate data collected by the lab so far, analyzing them and then making the data (and the visualizations) easily accessible to its users in a simple, click-to-operate user interface.

    ### The Data
    Data from the MWT Data Dashboard consolidated from a number of screens conducted by past students of the Rankin Lab. In these screens, students study straints of Caenorhabditis elegans carrying mutations (predominantly loss-of-function) in genes of interest.
    All of the experiments employ a 10-second ISI tap-habituation protocol. Statistics for the heatmap and the "phenomic profiles" of individual genes are done by comparing the strain of interest to the wildtype (N2) replicates conducted on the same day(s) of the experiment.

    ### Get Started
    This tool is designed to be easy to use - simply use your mouse to navigate to different elements of the dashboard to access the data. All the visualizations generated from the dashboard can be downloaded. 
    Navigate to appropraite pages to access relevant data. There are gene- and allele-specific pages if you are interested in looking at data for a specific gene (which is a mean of all the alleles analyzed) or allele.
    If you are interested in making comparisons between genes, navigate to the "Custome Gene/Allele Selection" pages.
    Not a fan of the artistic/design choices that were employed in these visualizations? Simple download the underlaying dataset and make your own graphs in your favourite language and visaulziation library.

    ### Contact Me
    Please reach me at joseph.liang@psych.ubc.ca for troubleshooting or feedback.
  

    """)
if page ==pages[6]:
    st.markdown("""
    ## References

    ### Multi-Worm Tracker:
    1. Swierczek, N. A., et al. (2011). High-throughput behavioral analysis in *C. elegans*. *Nature Methods, 8*(7), 592-598. https://doi.org/10.1038/nmeth.1625

    ### Neuronal Screen:
    2. Giles, A. C. (2012). *Candidate gene and high throughput genetic analysis of habituation in Caenorhabditis elegans*. University of British Columbia.

    ### G-Protein Screen:
    3. McEwan, A. (2013). *Modulation of habituation kinetics and behavioural shifts by members of the heterotrimeric G-protein signaling pathways*. University of British Columbia.

    ### Autism Spectrum Disorder (ASD) Screen:
    4. McDiarmid, T. A., et al. (2020). Systematic phenomics analysis of autism-associated genes reveals parallel networks underlying reversible impairments in habituation. *Proceedings of the National Academy of Sciences of the United States of America, 117*(1), 656-667. https://doi.org/10.1073/pnas.1912049116

    ### Unpublished Screens So Far:
    5. Liang, J. (n.d.). *Parkinson's disease screen*. Unpublished.

    6. Kepler, L. (n.d.). *Glia genes screen*. Unpublished.

    ### Online Resources:
    7. Harris, T. W., et al. (2020). WormBase: A modern model organism information resource. *Nucleic Acids Research, 48*(D1), D762-D767. https://doi.org/10.1093/nar/gkz920

    8. Kishore, R., et al. (2020). Automated generation of gene summaries at the Alliance of Genome Resources. *Database: The Journal of Biological Databases and Curation, 2020*, baaa037. https://doi.org/10.1093/database/baaa037

    ### CeNGEN:
    9. (Empty for now)
    """)

    
        

    
