# pages/gene.py
import streamlit as st
import io
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from utils.helpers import convert_df, read
from config import config

def render(data):
    st.header('Gene-specific Data')

    # Create a session state for the gene selection
    st.session_state.setdefault('gene_select', None)

    gene_list = [gene for gene in data["tap_output"]['Gene'].unique() if gene != 'N2']

    if st.session_state.gene_select in gene_list:
        default_index = gene_list.index(st.session_state.gene_select)
    else:
        default_index = 0 
        
    gene_option = st.selectbox(
        'Select a gene',
        gene_list,
        key="geneselect",
        index=default_index
    )
    st.session_state.gene_select = gene_option

    if gene_option:
        gene_id = data["id_data"].loc[data["id_data"]['Gene'] == gene_option, 'WBGene'].values
        if len(gene_id) == 0:
            gene_id = data["id_data"].loc[data["id_data"]['Sequence'] == gene_option, 'WBGene'].values

        if len(gene_id) > 0:
            glink = f'https://www.alliancegenome.org/gene/WB:{gene_id[0]}'
            st.markdown(f'<p style="font-size:20px">For more gene information on <a href="{glink}">{gene_option}</a> (Source: GenomeAlliance)</p>', unsafe_allow_html=True)
        else:
            st.markdown(f"information for {gene_option} not available")

    tap_output_gene = data["tap_output"][data["tap_output"]['Gene'] == gene_option]
    gene_tap_data = data["tap_output"][data["tap_output"]['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2', gene_option])].dropna(subset=['taps'])
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)

    col3, col4, col7 = st.columns([1, 1, 1])
    col3.subheader('Phenotypic profile')

    # seaborn plot
    sns.set_context('notebook', font_scale=1)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.barplot(x="Metric",
                     y="T_score", orient='v',
                     data=data["gene_profile_data"][data["gene_profile_data"].Gene == f"{gene_option}"],
                     palette=data["metric_palette"]).set_title(f"{gene_option}")
    plt.xticks(rotation=90)
    plt.ylabel("Normalized T-Score")
    plt.ylim(-3, 3)

    # download graph button
    gene_profile_plot = io.BytesIO()
    plt.savefig(gene_profile_plot, format='png', dpi=300, bbox_inches='tight')
    # display image
    col3.image(gene_profile_plot, width=None, caption=(f'Phenotypic profile of {gene_option}.'))
    col3_1, col3_2 = col3.columns(2)
    col3_1.download_button(label="Download Plot",
                        data=gene_profile_plot,
                        file_name=f"{gene_option}_profile.png",
                        mime="image/png",
                        key='dnldgeneprofile')
    col3_2.download_button(label="Download csv",
                        data=convert_df(data["gene_profile_data"][data["gene_profile_data"].Gene == f"{gene_option}"]),
                        file_name=f"Phenotypic profile of {gene_option}.csv",
                        mime="text/csv",
                        key='dnldgeneprofilecsv')

    col4.subheader('Rank in phenotype')
    gene_phenotype_option = col4.selectbox(
        'Select a phenotype',
        np.unique(data["phenotype_list"]),
        key='gene_phenotype_select')

    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    data_sorted = data["gene_MSD"].sort_values(by=[f"{gene_phenotype_option}-mean"])
    gene_colors = ["dimgray"] * len(data_sorted["Gene"])
    fig = go.Figure()

    # Scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            gene_colors = "red"
        elif row['Gene'] == gene_option:
            gene_colors = "magenta"
        else:
            gene_colors = "dimgray"
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
            showlegend=False,
            name=""
        ))

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
        yaxis=dict(showticklabels=True, dtick=1, tickfont=dict(color='black', size=6)),
        margin=dict(l=100, r=50, t=100, b=50),
        annotations=[
            dict(
                text=f'Sample mean distance from wildtype for all strains for selected phenotype: {gene_phenotype_option}. Error bars are 95% CI',
                xref="paper", yref="paper",
                x=0, y=-0.2,
                showarrow=False,
                font=dict(size=12, color="black")
            )
        ]
    )


    # recreate matplotlib chart for download button 
    fig_mpl, ax = plt.subplots(figsize=(6, 12))

    # Loop over each row to plot individual points with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            color = "red"
        elif row['Gene'] == gene_option:
            color = "magenta"
        else:
            color = "dimgray"

        ax.errorbar(
            x=row[f"{gene_phenotype_option}-mean"],
            y=row["Gene"],
            xerr=[[row[f"{gene_phenotype_option}-mean"] - row[f"{gene_phenotype_option}-ci95_lo"]],
                [row[f"{gene_phenotype_option}-ci95_hi"] - row[f"{gene_phenotype_option}-mean"]]],
            fmt='o',
            color=color,
            ecolor=color,
            elinewidth=1,
            capsize=3,
            markersize=6
        )

    # Vertical reference line and styling
    ax.axvline(x=0, color='red', linestyle='--')
    ax.set_xlabel("Sample Mean Distance")
    ax.set_ylabel("Gene")
    ax.set_title(f"{gene_phenotype_option}", fontsize=14)
    ax.set_yticks(range(len(data_sorted["Gene"])))
    ax.set_yticklabels(data_sorted["Gene"], fontsize=6)
    ax.invert_yaxis()
    plt.tight_layout()

    # Add annotation below chart
    plt.figtext(
        0.5, -0.05,
        f'Sample mean distance from wildtype for all strains for selected phenotype: {gene_phenotype_option}. Error bars are 95% CI',
        wrap=True, ha='center', fontsize=10
    )

    # Save to buffer for download
    gene_phenotype_plot = io.BytesIO()
    plt.savefig(gene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    gene_phenotype_plot.seek(0)
    plt.close()


    col4.plotly_chart(fig, use_container_width=True, **{'config': config})

    # UPDATED: using combined data
    gene_dat = data_sorted[["Gene", f"{gene_phenotype_option}-mean", f"{gene_phenotype_option}-ci95_lo", f"{gene_phenotype_option}-ci95_hi"]]
    gene_dat.columns = ["Gene", f"{gene_phenotype_option}", f"{gene_phenotype_option}-lower", f"{gene_phenotype_option}-upper"]

    # Insert download graph button
    col4_1, col4_2 = col4.columns(2)
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
                    data=gene_tap_data_plot,
                    hue='Gene', # Here we use the extra column from step 6 to separate by group
                    palette=["black" if gene == "N2" else "darkorange" for gene in gene_tap_data_plot['Gene'].unique()],
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
                tabs[i].image(img1, width=None, caption=f'Habituation of Response {metric}: {gene_option}')
                # Insert download plot and download csv button
                col1, col2 = tabs[i].columns(2)
                col1.download_button("Download Plot", data=img1, file_name=f"{metric} of Tap Habituation {gene_option}.png", mime="image/png", key=f'dnldbtngene_{i}')
                col2.download_button("Download csv", data=convert_df(gene_tap_data_plot), file_name=f"Gene-specific Data {gene_option}.csv", mime="text/csv", key=f'dnldbtngene2_{i}')


    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data', key='readgenebaseoutcsv'):
        read_data_flag = True
        st.warning("Please wait, another button will open up to download the data. This might take several minutes ")

    # If the button is pressed, read the data and then show show button to download it
    if read_data_flag:
        baseline_output = read('tap_baseline_data')
        baseline_output = baseline_output[baseline_output['Screen'].isin(data["datasets"])].replace(["N2_N2", "N2_XJ1"], "N2")
        st.download_button(
            label="Download raw baseline data",
            data=convert_df(baseline_output),
            file_name=f"raw_baseline_data.csv",
            mime="text/csv",
            key='dnldgenebaseoutcsv')
