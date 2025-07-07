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
metric_palette = [
    "k", "k", "k",
    "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgrey", "darkgray",
    "lightsteelblue", "lightsteelblue", "lightsteelblue",
    "powderblue", "powderblue", "powderblue",
    "cadetblue", "cadetblue", "cadetblue",
    "thistle", "thistle", "thistle",
    "slateblue", "slateblue", "slateblue"
]
