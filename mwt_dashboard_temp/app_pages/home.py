# pages/home.py
import streamlit as st
import io
import plotly.graph_objects as go
import numpy as np
from utils.helpers import convert_df, read
from config import config

def render(data):
    st.header('Home - Data at a Glance')

    # Visualisations for data tab
    col1, col2 = st.columns([4, 5])

    col1.subheader("For A Single Phenotype")

    phenotype_option = col1.selectbox(
        'Select a phenotype',
        np.unique(data["phenotype_list"]), key="phenotypeselect")

    data_sorted = data["gene_MSD"].sort_values(by=[f"{phenotype_option}-mean"]).reset_index(drop=True)
    
    colors = ["dimgrey"] * len(data_sorted)

    fig = go.Figure()

    # Add scatter plot with error bars
    for i, row in data_sorted.iterrows():
        if row['Gene'] == "N2":
            colors = "red"
        else:
            colors = "dimgray"
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
            showlegend=False,
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
        margin=dict(l=10, r=10, t=40, b=0),
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
    data_dat = data["gene_MSD"].sort_values(by=[f"{phenotype_option}-mean"])[["Gene", f"{phenotype_option}-mean", f"{phenotype_option}-ci95_lo", f"{phenotype_option}-ci95_hi"]]
    data_dat.columns = ["Gene", f"{phenotype_option}", f"{phenotype_option}-lower", f"{phenotype_option}-upper"]

    # Insert download graph button
    col1_1, col1_2 = col1.columns(2)
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
        z=data["tap_tstat_data"].set_index('Gene').values,
        x=data["tap_tstat_data"].set_index('Gene').columns,
        y=data["tap_tstat_data"].set_index('Gene').index,
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

    col2_1, col2_2 = col2.columns(2)
    col2_1.download_button(label="Download Plot",
                        data=imgheatmap,
                        file_name="Heatmap.png",
                        mime="image/png",
                        key='dnldheatmap')
    # Add download buttons    
    col2_2.download_button(
        label="Download CSV",
        data=convert_df(data["tap_tstat_data"].set_index('Gene')),
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
        baseline_output = baseline_output[baseline_output['Screen'].isin(data["datasets"])].replace(["N2_N2", "N2_XJ1"], "N2")

        st.download_button(label="Download raw baseline data",
                        data=convert_df(baseline_output),
                        file_name=f"raw_baseline_data.csv",
                        mime="text/csv",
                        key='dnldbaseoutcsv')
