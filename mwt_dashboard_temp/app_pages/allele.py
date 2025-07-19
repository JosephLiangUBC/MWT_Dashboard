# pages/allele.py
import streamlit as st
import io
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from utils.helpers import convert_df, read
from config import config

def render(data):
    st.header("Allele-specific Data")

    # Create a session state for the allele selection
    st.session_state.setdefault('allele_select', None)

    allele_list = [allele for allele in data["tap_output"]['dataset'].unique() if allele != 'N2']

    if st.session_state.allele_select in allele_list:
        default_index = allele_list.index(st.session_state.allele_select)
    else:
        default_index = 0 

    allele_option = st.selectbox(
        'Select an allele',
        allele_list, 
        key="alleleselect",
        index=default_index
    )
    st.session_state.allele_select = allele_option

    # splititing gene allele to get gene and allele columns
    gene, allele = allele_option.split('_')

    # check if allele and gene options match the preexisting list
    if allele_option:
        gene_id = data["id_data"].loc[data["id_data"]['Gene'] == gene, 'WBGene']
        allele_id = data["id_data"].loc[
            (data["id_data"]['Gene'] == gene) &
            (data["id_data"]['Allele'] == allele_option),
            'WBAllele'
        ]

        allele_id = allele_id.values[0] if allele_id.any() else "default_value"
        gene_id = gene_id.values[0] if gene_id.any() else "default_value"

        # links to the Alliance Genome and WormBase website 
        glink = f'https://www.alliancegenome.org/gene/WB:{gene_id}'
        wlink = f'https://wormbase.org/species/c_elegans/variation/{allele_id}'

    # display links 
    st.markdown(
        f'<p style="font-size:20px">For more gene information on <a href="{glink}">{gene}</a> (Source: GenomeAlliance) and allele information on <a href="{wlink}">{allele_option}</a> (Source: WormBase)</p>',
        unsafe_allow_html=True
        )

    allele_profile_data = data["allele_profile_data"]
    phenotype_list = data["phenotype_list"]

    tap_output_allele = data["tap_output"][data["tap_output"]['dataset'] == allele_option]
    allele_tap_data = data["tap_output"][data["tap_output"]['Date'].isin(tap_output_allele['Date'].unique())]
    allele_tap_data_plot = allele_tap_data[allele_tap_data['dataset'].isin(['N2', allele_option])].dropna(subset=['taps'])
    allele_tap_data_plot['taps'] = allele_tap_data_plot['taps'].astype(int)

    col3, col4, col7 = st.columns([1, 1, 1])
    col3.subheader('Phenotypic profile')

    # seaborn plot
    sns.set_context('notebook', font_scale=1)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.barplot(x="Metric", # <- Here we use seaborn as our graphing package.
                        y="T_score", orient='v',
                        data=allele_profile_data[allele_profile_data.dataset == f"{allele_option}"],
                        palette=data["metric_palette"]).set_title(f"{allele_option}")
    plt.xticks(rotation=90)
    plt.ylabel("Normalized T-Score") # <- X-axis title
    plt.ylim(-3, 3)

    # download graph button
    allele_profile_plot = io.BytesIO()
    plt.savefig(allele_profile_plot, format='png', dpi=300, bbox_inches='tight')
    # display image 
    col3.image(allele_profile_plot, width=None, caption=(f'Phenotypic profile of {allele_option}.'))
    # download button for plots
    col3_1, col3_2 = col3.columns(2)
    col3_1.download_button(label="Download Plot",
                            data=allele_profile_plot,
                            file_name=f"{allele_option}_profile.png",
                            mime="image/png",
                            key='dnldalleleprofile')
    # download button for data
    col3_2.download_button(label="Download csv",
                            data=convert_df(allele_profile_data[allele_profile_data.dataset == f"{allele_option}"]),
                            file_name=f"Phenotypic profile of gene-allele {allele_option}.csv",
                            mime="text/csv",
                            key='dnldalleleprofilecsv')

    col4.subheader('Rank in phenotype')

    allele_phenotype_option = col4.selectbox(
        'Select a phenotype',
        np.unique(phenotype_list),
        key='allele_phenotype_select'
    )
    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot

    data_sorted = data["allele_MSD"].sort_values(by=[f"{allele_phenotype_option}-mean"])
    allele_colors = ["dimgray"] * len(data_sorted["dataset"])
    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['dataset'] == "N2":
            allele_colors = "red"
        elif row['dataset'] == allele_option:
            allele_colors = "magenta"
        else:
            allele_colors = "dimgray"
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
            showlegend=False, # Hide individual points from legend
            name=""
        ))

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
        yaxis=dict(showticklabels=True, dtick=1, tickfont=dict(color='black', size=6)),
        margin=dict(l=100, r=50, t=100, b=50), # Adjust margins
        annotations=[
            dict(
                text=f'Sample mean distance from wildtype for all strains for selected phenotype: {allele_phenotype_option}. Error bars are 95% CI',
                xref="paper", yref="paper",
                x=0, y=-0.2,
                showarrow=False,
                font=dict(size=12, color="black")
            )
        ]
    )

    allele_phenotype_plot = io.BytesIO()
    fig.write_image(allele_phenotype_plot, format='png', scale=3)
    allele_phenotype_plot.seek(0)
    col4.plotly_chart(fig, use_container_width=True, **{'config': config})

    allele_dat = data_sorted[["dataset", f"{allele_phenotype_option}-mean", f"{allele_phenotype_option}-ci95_lo", f"{allele_phenotype_option}-ci95_hi"]]
    allele_dat.columns = ["gene-allele", f"{allele_phenotype_option}", f"{allele_phenotype_option}-lower", f"{allele_phenotype_option}-upper"]

    # Insert download graph button
    col4_1, col4_2 = col4.columns(2)
    col4_1.download_button(label="Download Plot",
                            data=allele_phenotype_plot,
                            file_name=f"{allele_option}_{allele_phenotype_option}_profile.png",
                            mime="image/png",
                            key='dnldallelephenotypeprofile')
    col4_2.download_button(label="Download csv",
                            data=convert_df(allele_dat),
                            file_name=f"Allele-specific Data Sample mean distance {allele_phenotype_option}.csv",
                            mime="text/csv",
                            key='dnldallelephenotypeprofilecsv')

    # Habituation Curves 
    col7.subheader('Habituation Curves of Response')

    with col7:

        metrics = [ "Probability", "Duration", "Speed",
                   "PSA Instantaneous Speed", "PSA Interval Speed",
                    "PSA Bias", "PSA Kink", "PSA Crab",
                    "PSA Aspect Ratio", "PSA Curve"]
        
        tabs = st.tabs(metrics)

        for i, metric in enumerate(metrics):
            with tabs[i]:
                #  Habituation Plot
                fig, ax = plt.subplots(figsize=(12, 10))
                plt.gca().xaxis.grid(False) # gets rid of x-axis markers to make data look clean
                sns.pointplot(
                    x="taps", # Here we use seaborn as our graphing package.
                    y=metric,
                    data=allele_tap_data_plot,
                    hue='dataset', # Here we use the extra column from step 6 to separate by group
                    palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
                    errorbar='se' # Confidence interval. 95 = standard error
                )
                plt.xlabel("Taps", fontsize='12') # x-axis title
                plt.ylabel(metric, fontsize='12') # y-axis title
                plt.title(f"Habituation of Response {metric}", fontsize='16') # figure title
                plt.ylim(0, 1 if metric == "Probability" else None)
                ax.legend(loc='upper right', fontsize='12') # legend location

                # download graph button
                img1 = io.BytesIO()
                plt.savefig(img1, format='png', dpi=300, bbox_inches='tight')
                img1.seek(0)
                # display image
                tabs[i].image(img1, width=None, caption=f'Habituation of Response {metric}: {allele_option}')
                # Insert download plot and download csv button
                col1, col2 = tabs[i].columns(2)
                col1.download_button("Download Plot", data=img1, file_name=f"{metric} of Tap Habituation {allele_option}.png", mime="image/png", key=f'dnldbtngene_{i}')
                col2.download_button("Download csv", data=convert_df(allele_tap_data_plot), file_name=f"Gene-specific Data {allele_option}.csv", mime="text/csv", key=f'dnldbtngene2_{i}')

        







    # with col7:
    #     tab1, tab2, tab3 = st.tabs(["Probability", "Duration", "Speed"])

    #     with tab1:
    #         fig, ax = plt.subplots(figsize=(12, 10))
    #         plt.gca().xaxis.grid(False) # gets rid of x-axis markers to make data look clean
    #         ax = sns.pointplot(x="taps", # Here we use seaborn as our graphing package.
    #                             y="prob",
    #                             data=allele_tap_data_plot,
    #                             hue='dataset', # Here we use the extra column from step 6 to separate by group
    #                             palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
    #                             errorbar='se') # Confidence interval. 95 = standard error
    #         plt.xlabel("Taps")
    #         plt.ylabel("Probability")
    #         plt.title("Habituation of Response Probability", fontsize='16')
    #         plt.ylim(0, 1)
    #         ax.legend(loc='upper right', fontsize='12')

    #         img1 = io.BytesIO()
    #         plt.savefig(img1, format='png', dpi=300, bbox_inches='tight')
    #         tab1.image(img1, width=None, caption=f'Habituation of Response Probability: {allele_option}')
    #         # Insert download plot and download csv button
    #         tab1_1, tab1_2 = tab1.columns(2)
    #         tab1_1.download_button(label="Download Plot",
    #                                 data=img1,
    #                                 file_name=f"Probability of Tap Habituation {allele_option}.png",
    #                                 mime="image/png",
    #                                 key='dnldbtn1')
    #         tab1_2.download_button(label="Download csv",
    #                                 data=convert_df(allele_tap_data_plot),
    #                                 file_name=f"Allele-specific Data {allele_option}.csv",
    #                                 mime="text/csv",
    #                                 key='dnldbtn7')

    #     with tab2:
    #         #  Habituation of Response Duration Plot
    #         fig, ax = plt.subplots(figsize=(12, 10))
    #         # seaborn plot
    #         ax = sns.pointplot(x="taps",
    #                             y="dura",
    #                             data=allele_tap_data_plot,
    #                             hue='dataset',
    #                             palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
    #                             errorbar='se')
    #         plt.xlabel("Taps", fontsize='12')
    #         plt.ylabel("Duration", fontsize='12')
    #         plt.title("Habituation of Response Duration", fontsize='16')
    #         plt.ylim(0, None)
    #         ax.legend(loc='upper right', fontsize='12')

    #         img2 = io.BytesIO()
    #         plt.savefig(img2, format='png', dpi=300, bbox_inches='tight')
    #         tab2.image(img2, width=None, caption=f'Habituation of Response Duration: {allele_option}')
    #         # Insert download plot and download csv button
    #         tab2_1, tab2_2 = tab2.columns(2)
    #         tab2_1.download_button(label="Download Plot",
    #                                 data=img2,
    #                                 file_name=f"Duration of Tap Habituation {allele_option}.png",
    #                                 mime="image/png",
    #                                 key='dnldbtn5')
    #         tab2_2.download_button(label="Download csv",
    #                                 data=convert_df(allele_tap_data_plot),
    #                                 file_name=f"Allele-specific Data {allele_option}.csv",
    #                                 mime="text/csv",
    #                                 key='dnldbtn11')

    #     with tab3:
    #         #  Habituation of Response Speed Plot
    #         fig, ax = plt.subplots(figsize=(12, 10))
    #         # seaborn plot
    #         ax = sns.pointplot(x="taps",
    #                             y="speed",
    #                             data=allele_tap_data_plot,
    #                             hue='dataset',
    #                             palette=["black" if gene == "N2" else "darkorange" for gene in allele_tap_data_plot['dataset'].unique()],
    #                             errorbar='se')
    #         plt.xlabel("Taps", fontsize='12')
    #         plt.ylabel("Speed", fontsize='12')
    #         plt.title("Habituation of Response Speed", fontsize='16')
    #         plt.ylim(0, None)
    #         ax.legend(loc='upper right', fontsize='12')

    #         img3 = io.BytesIO()
    #         plt.savefig(img3, format='png', dpi=300, bbox_inches='tight')
    #         tab3.image(img3, width=None, caption=f'Habituation of Response Speed: {allele_option}')
    #         tab3_1, tab3_2 = tab3.columns(2)
    #         tab3_1.download_button(label="Download Plot",
    #                                 data=img3,
    #                                 file_name=f"Speed of Tap Habituation {allele_option}.png",
    #                                 mime="image/png",
    #                                 key='dnldbtn6')
    #         tab3_2.download_button(label="Download csv",
    #                                 data=convert_df(allele_tap_data_plot),
    #                                 file_name=f"Allele-specific Data {allele_option}.csv",
    #                                 mime="text/csv",
    #                                 key='dnldbtn12')



    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data',key='readallelebaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(data["datasets"])].replace(["N2_N2", "N2_XJ1"], "N2")
        conn.close()
        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldallelebaseoutcsv')
