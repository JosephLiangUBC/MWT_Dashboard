# config.py

# This controls the page names shown in the sidebar
pages = [
    "Home - Data at a Glance",
    "Gene-specific Data",
    "Allele-specific Data",
    "Custom Gene Selection",
    "Custom Allele Selection",
    "Post Stimulus Data", 
    "Help & Documentation",
    "Citations"
]

# Plotly config for all download buttons
config = {
    'toImageButtonOptions': {
        'format': 'png',
        'filename': 'plot',
        'height': None,
        'width': None,
        'scale': 3
    }
}

# color palette used in barplots
# metric_palette = [
#     "k", "k", "k",
#     "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgrey", "darkgray",
#     "lightsteelblue", "lightsteelblue", "lightsteelblue",
#     "powderblue", "powderblue", "powderblue",
#     "cadetblue", "cadetblue", "cadetblue",
#     "thistle", "thistle", "thistle",
#     "slateblue", "slateblue", "slateblue"
# ]

metric_palette=["k","k","k",
                "darkgray","darkgray","darkgray","darkgray","darkgray","darkgray","darkgray","darkgray",
                "lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue","lightsteelblue",
                "powderblue","powderblue","powderblue","powderblue","powderblue","powderblue","powderblue","powderblue","powderblue","powderblue",
                "orange","orange","orange","orange","orange","orange","orange",
                "darkorange","darkorange","darkorange","darkorange","darkorange","darkorange","darkorange",
                "goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod","goldenrod",
                "cadetblue","cadetblue","cadetblue","cadetblue","cadetblue","cadetblue","cadetblue","cadetblue","cadetblue","cadetblue",
                "thistle","thistle","thistle","thistle","thistle","thistle","thistle","thistle","thistle","thistle",
                "slateblue","slateblue","slateblue","slateblue","slateblue","slateblue","slateblue","slateblue","slateblue","slateblue",
                "lightcoral","lightcoral","lightcoral","lightcoral","lightcoral","lightcoral","lightcoral",
                ]