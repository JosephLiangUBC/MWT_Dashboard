# pages/custom_allele.py
import streamlit as st
import io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import itertools
import sqlite3
from utils.helpers import convert_df, read

def render(data):
    st.header('Custom Allele Selection')
    st.session_state.setdefault('allele_select', [allele for allele in data["tap_output"]['dataset'].unique() if allele != 'N2'][0])

    allele_multiple = st.multiselect(
        label="Select Allele",
        options=[allele for allele in data["tap_output"]['dataset'].unique() if allele != 'N2'],
        default=st.session_state.allele_select,
        placeholder="make a selection",
        help="select and de-select alleles you want to analyze",
        key="alleleselection")
    st.session_state.allele_select = allele_multiple

    na_list = []
    g_link_list = [] # list for gene links (AllianceGenome)
    w_link_list = [] # list for allele links (WormBase)
    allele_list = []
    for a in allele_multiple:
        gene, allele = a.split('_')
        allele_list.append(a)
        gene_id = data["id_data"].loc[data["id_data"]['Gene'] == gene, 'WBGene']
        allele_id = data["id_data"].loc[(data["id_data"]['Gene'] == gene) & (data["id_data"]['Allele'] == allele), 'WBAllele']

        if allele_id.any():
            allele_id = allele_id.values[0]
            wlink = f'https://wormbase.org/species/c_elegans/variation/{allele_id}'
            w_link_list.append(f'<a href="{wlink}">{allele}</a>')
        else: # if allele option doesnt match then send to default page
            na_list.append(allele)

        if gene_id.any():
            gene_id = gene_id.values[0]
            glink = f'https://www.alliancegenome.org/gene/WB:{gene_id}'
            g_link_list.append(f'<a href="{glink}">{gene}</a>')
        else:
            gene_id= "default value"
            na_list.append(gene)

    st.markdown(f"<p style='font-size:20px'>For more gene information on {', '.join(g_link_list)} (Source: GenomeAlliance) </p>", unsafe_allow_html=True)
    st.markdown(f"<p style='font-size:20px'>For more allele information on {', '.join(w_link_list)} (Source: WormBase) </p>", unsafe_allow_html=True)
    if na_list:
        na_links = [f'<a href="https://wormbase.org/species/c_elegans/variation/">{allele}</a>' for allele in na_list]
        st.markdown(f"<p style='font-size:20px'>Information not available for: {', '.join(na_links)}</p>", unsafe_allow_html=True)

    #filter data for particular allele
    tap_output_allele = data["tap_output"][data["tap_output"]['dataset'].isin(allele_multiple)]
    allele_tap_data = data["tap_output"][data["tap_output"]['Date'].isin(tap_output_allele['Date'].unique())]
    allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2'] + allele_multiple)].dropna(subset=['taps'])
    allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)

    #add columns for msd, habituation plots and heatmap plots
    col12, col13, col14 = st.columns([1, 1, 1])
    # Filter the dataframe for the selected genes
    tap_tstat_allele_selected = data["tap_tstat_allele"][data["tap_tstat_allele"]['dataset'].isin(allele_list)]

    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_allele_selected.set_index('dataset').values,
        x=tap_tstat_allele_selected.set_index('dataset').columns,
        y=tap_tstat_allele_selected.set_index('dataset').index,
        colorscale='RdBu',
        zmin=-3,
        zmax=3,
        colorbar=dict(
            len=0.95,
            thickness=10,
            tickvals=[-3, 0, 3],
            ticktext=['-3', '0', '3'],
            title="",
            # titleside="right"
        )
    ))
    h = 200 * len(allele_multiple)
    fig.update_layout(
        width=900,
        height=h,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(tickangle=0),
        xaxis=dict(showticklabels=True, tickfont=dict(size=8))
    )

    imgheatmap = io.BytesIO()
    fig.write_image(imgheatmap, format='png', scale=3)
    imgheatmap.seek(0)

    col12.subheader(f'Comprehensive heatmap of the dataset with selected alleles')
    col12.plotly_chart(fig, use_container_width=True, config=data["plotly_config"])

    col12_1, col12_2 = col12.columns(2)
    col12_1.download_button(label="Download Plot",
                            data=imgheatmap,
                            file_name="Heatmap.png",
                            mime="image/png",
                            key='dnldheatmapcustomallele')
    col12_2.download_button(label="Download CSV",
                            data=convert_df(tap_tstat_allele_selected.set_index('dataset')),
                            file_name="Data_Glance_Heatmap.csv",
                            mime="text/csv",
                            key='dnldheatmapcsvcustomallele')

    col13.subheader('Rank in phenotype')
    multiallele_phenotype_option = col13.selectbox(
        'Select a phenotype',
        np.unique(data["phenotype_list"]),
        key='multiallele_phenotype_select')
    data_sorted = data["allele_MSD"].sort_values(by=[f"{multiallele_phenotype_option}-mean"])
    allele_colors = ["dimgray"] * len(data_sorted["dataset"])
    ticktext = []
    tickvals = []
    fig = go.Figure()
    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['dataset'] == "N2":
            allele_colors = "red"
            ticktext.append(row['dataset'])
            tickvals.append(row['dataset'])
        elif row['dataset'] in allele_multiple:
            allele_colors = "magenta"
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
        margin=dict(l=100, r=50, t=100, b=50),
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
    fig.write_image(multiallele_phenotype_plot, format='png', scale=3)
    multiallele_phenotype_plot.seek(0)
    col13.plotly_chart(fig, use_container_width=True, config=data["plotly_config"])

    #combine data and rename columns :
    multiallele_dat = pd.concat([
        data["allele_MSD"].sort_values(by=[f"{multiallele_phenotype_option}-mean"])["dataset"],
        data["allele_MSD"].sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-mean"],
        data["allele_MSD"].sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-ci95_lo"],
        data["allele_MSD"].sort_values(by=[f"{multiallele_phenotype_option}-mean"])[f"{multiallele_phenotype_option}-ci95_hi"]],
        axis=1)
    multiallele_dat.columns = ["Allele", f"{multiallele_phenotype_option}", f"{multiallele_phenotype_option}-lower", f"{multiallele_phenotype_option}-upper"]

    col13_1, col13_2 = col13.columns(2)
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
    colors_list = sns.color_palette("husl", n_colors=len(alleles) + 1)
    color_cycle = itertools.cycle(colors_list)
    colors = [next(color_cycle) for _ in range(len(alleles))]
    new_palette = ["black" if allele == "N2" else color for allele, color in zip(alleles, colors)]

    with col14:
        tab10, tab11, tab12 = st.tabs(["Probability", "Duration", "Speed"])
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
    if st.button('Read and get Download button for Baseline Data', key='readallelemultibaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        for conn_path in data["conn_list"]:
            try:
                conn = sqlite3.connect(conn_path)
                break
            except:
                pass
        baseline_output = read('tap_baseline_data', conn)
        baseline_output = baseline_output[baseline_output['Screen'].isin(data["datasets"])].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                           data=convert_df(baseline_output[baseline_output['dataset'].isin(allele_multiple)]),
                           file_name=f"raw_baseline_data.csv",
                           mime="text/csv",
                           key='dnldallelemultibaseoutcsv')