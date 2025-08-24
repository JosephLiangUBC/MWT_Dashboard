# pages/custom_gene.py
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
    # multiple selection option for genes
    st.header('Custom Gene Selection')
    st.session_state.setdefault('gene_select', [gene for gene in data["tap_output"]['Gene'].unique() if gene != 'N2'][0])

    gene_multiple = st.multiselect(
        label="Select Genes",
        options=[gene for gene in data["tap_output"]['Gene'].unique() if gene != 'N2'],
        default=st.session_state.gene_select,
        placeholder="make a selection",
        help="select and de-select genes you want to analyze",
        key="geneselection")
    st.session_state.gene_select = gene_multiple

    na_list = []
    g_link_list = []
    for gene in gene_multiple:
        gene_id = data["id_data"].loc[data["id_data"]['Gene'] == gene, 'WBGene'].values
        if len(gene_id) == 0:
            gene_id = data["id_data"].loc[data["id_data"]['Sequence'] == gene, 'WBGene'].values
        if len(gene_id) > 0:
            glink = f'https://www.alliancegenome.org/gene/WB:{gene_id[0]}'
            g_link_list.append(f'<a href="{glink}">{gene}</a>')
        else:
            na_list.append(gene)
    st.markdown(f"<p style='font-size:20px'>For more gene information on {', '.join(g_link_list)} (Source: GenomeAlliance)</p>", unsafe_allow_html=True)
    if na_list:
        na_links = [f'<a href="https://www.alliancegenome.org">{gene}</a>' for gene in na_list]
        st.markdown(f"<p style='font-size:20px'>Information not available for: {', '.join(na_links)}</p>", unsafe_allow_html=True)

    # filter data for particular genes
    tap_output_gene = data["tap_output"][data["tap_output"]['Gene'].isin(gene_multiple)]
    gene_tap_data = data["tap_output"][data["tap_output"]['Date'].isin(tap_output_gene['Date'].unique())]
    gene_tap_data_plot = gene_tap_data[gene_tap_data['Gene'].isin(['N2'] + gene_multiple)].dropna(subset=['taps'])
    gene_tap_data_plot['taps'] = gene_tap_data_plot['taps'].astype(int)

    col9, col10, col11 = st.columns([1, 1, 1])
    #current
    tap_tstat_selected = data["tap_tstat_data"][data["tap_tstat_data"]['Gene'].isin(gene_multiple)]

    # Create a heatmap
    fig = go.Figure(data=go.Heatmap(
        z=tap_tstat_selected.set_index('Gene').values,
        x=tap_tstat_selected.set_index('Gene').columns,
        y=tap_tstat_selected.set_index('Gene').index,
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
    h = 200 * len(gene_multiple)
    fig.update_layout(
        width=900,
        height=h,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_title="",
        yaxis_title="",
        yaxis=dict(tickangle=0),
        xaxis=dict(showticklabels=True, tickfont=dict(size=8))
    )

    fig_mpl, ax = plt.subplots(figsize=(9, 0.2 * len(gene_multiple)))

    heatmap_data = tap_tstat_selected.set_index('Gene')
    cax = ax.imshow(heatmap_data.values, cmap='RdBu', aspect='auto', vmin=-3, vmax=3)

    # Tick labels
    ax.set_xticks(range(len(heatmap_data.columns)))
    ax.set_xticklabels(heatmap_data.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(heatmap_data.index)))
    ax.set_yticklabels(heatmap_data.index, fontsize=8)

    # Colorbar
    cbar = fig_mpl.colorbar(cax, ax=ax, orientation='vertical', fraction=0.025, pad=0.02)
    cbar.set_ticks([-3, 0, 3])
    cbar.set_ticklabels(['-3', '0', '3'])

    plt.tight_layout()

    # Save to BytesIO
    imgheatmap = io.BytesIO()
    plt.savefig(imgheatmap, format='png', dpi=300, bbox_inches='tight')
    imgheatmap.seek(0)

    col9.subheader('Comprehensive heatmap of the dataset with selected genes')
    col9.plotly_chart(fig, use_container_width=True, config=data["plotly_config"])

    col9_1, col9_2 = col9.columns(2)
    col9_1.download_button(label="Download Plot",
                            data=imgheatmap,
                            file_name="Heatmap.png",
                            mime="image/png",
                            key='dnldheatmapcustom')
    col9_2.download_button(label="Download CSV",
                            data=convert_df(tap_tstat_selected.set_index('Gene')),
                            file_name="Data_Glance_Heatmap.csv",
                            mime="text/csv",
                            key='dnldheatmapcsvcustom')

    col10.subheader('Rank in phenotype')
    multigene_phenotype_option = col10.selectbox(
        'Select a phenotype',
        np.unique(data["phenotype_list"]),
        key='multigene_phenotype_select')
    # seaborn graph of phenotypic view (sample mean distance) + st.pyplot
    data_sorted = data["gene_MSD"].sort_values(by=[f"{multigene_phenotype_option}-mean"])
    gene_colors = ["dimgray"] * len(data_sorted["Gene"])
    ticktext = []
    tickvals = []
    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            gene_colors = "red"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        elif row['Gene'] in gene_multiple:
            gene_colors = "magenta"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        else:
            gene_colors = "dimgray"
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
            showlegend=False, # Hide individual points from legend
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
        margin=dict(l=100, r=50, t=100, b=50),
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
    
    # recreate matplotlib chart for download button
    fig_mpl, ax = plt.subplots(figsize=(6, 12))
    ticktext = []
    tickvals = []
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            gene_colors = "red"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        elif row['Gene'] in gene_multiple:
            gene_colors = "magenta"
            ticktext.append(row['Gene'])
            tickvals.append(row['Gene'])
        else:
            gene_colors = "dimgray"

        ax.errorbar(
            x=row[f"{multigene_phenotype_option}-mean"],
            y=row["Gene"],
            xerr=[[row[f"{multigene_phenotype_option}-mean"] - row[f"{multigene_phenotype_option}-ci95_lo"]],
                [row[f"{multigene_phenotype_option}-ci95_hi"] - row[f"{multigene_phenotype_option}-mean"]]],
            fmt='o',
            color=gene_colors,
            ecolor=gene_colors,
            elinewidth=1,
            capsize=3
        )

    ax.axvline(x=0, color='red', linestyle='--')
    ax.set_title(f"{multigene_phenotype_option}", fontsize=14)
    ax.set_xlabel("Sample Mean Distance")
    ax.set_ylabel("Gene")
    ax.set_yticks(tickvals)
    ax.set_yticklabels(ticktext)
    ax.invert_yaxis()
    plt.tight_layout()

    plt.figtext(
        0.5, -0.05,
        f'Sample mean distance from wildtype for all strains for selected phenotypes: {multigene_phenotype_option}. Error bars are 95% CI',
        wrap=True, ha='center', fontsize=10
    )

    multigene_phenotype_plot = io.BytesIO()
    plt.savefig(multigene_phenotype_plot, format='png', dpi=300, bbox_inches='tight')
    multigene_phenotype_plot.seek(0)
    plt.close()

    col10.plotly_chart(fig, use_container_width=True, config=data["plotly_config"])

    #combine data and rename columns :
    multigene_dat = pd.concat([
        data["gene_MSD"].sort_values(by=[f"{multigene_phenotype_option}-mean"])["Gene"],
        data["gene_MSD"].sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-mean"],
        data["gene_MSD"].sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-ci95_lo"],
        data["gene_MSD"].sort_values(by=[f"{multigene_phenotype_option}-mean"])[f"{multigene_phenotype_option}-ci95_hi"]],
        axis=1)
    multigene_dat.columns = ["Gene", f"{multigene_phenotype_option}", f"{multigene_phenotype_option}-lower", f"{multigene_phenotype_option}-upper"]

    # Insert download graph button
    col10_1, col10_2 = col10.columns(2)
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
    colors_list = sns.color_palette("husl", n_colors=len(genes) + 1)
    color_cycle = itertools.cycle(colors_list)

    # Create a list of unique colors for the alleles
    colors = [next(color_cycle) for _ in range(len(genes))]

    # Create a palette with 'black' for 'N2' and the unique colors for the other genes
    new_palette = ["black" if gene == "N2" else color for gene, color in zip(genes, colors)]

    with col11:

        metrics = [ "Probability", "Duration", "Speed",
                   "PSA Speed", 
                #    "PSA Interval Speed",
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
                    palette=new_palette,
                    errorbar='se' # Confidence interval. 95 = standard error
                )
                plt.xlabel("Taps", fontsize='12') # x-axis title
                plt.ylabel(metric, fontsize='12') # y-axis title
                plt.title(f"Habituation of Response {metric}", fontsize='16') # figure title
                plt.ylim(0, 1 if metric == "Probability" else None)
                ax.legend(loc='upper right', fontsize='12') # legend location

                # download graph button
                img7 = io.BytesIO()
                plt.savefig(img7, format='png', dpi=300, bbox_inches='tight')
                img7.seek(0)
                # display image
                tabs[i].image(img7, width=None, caption=f'Habituation of Response {metric}: {gene_multiple}')
                # Insert download plot and download csv button
                col1, col2 = tabs[i].columns(2)
                col1.download_button("Download Plot", data=img7, file_name=f"{metric} of Tap Habituation {gene_multiple}.png", mime="image/png", key=f'dnldbtncustgene_{i}')
                col2.download_button("Download csv", data=convert_df(gene_tap_data_plot), file_name=f"Gene-specific Data {gene_multiple}.csv", mime="text/csv", key=f'dnldbtncustgene2_{i}')


    # Create a flag variable
    read_data_flag = False
    if st.button('Read and get Download button for Baseline Data', key='readgenemultibaseoutcsv'):
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
                           data=convert_df(baseline_output[baseline_output['Gene'].isin(gene_multiple)]),
                           file_name=f"raw_baseline_data.csv",
                           mime="text/csv",
                           key='dnldgenemultibaseoutcsv')