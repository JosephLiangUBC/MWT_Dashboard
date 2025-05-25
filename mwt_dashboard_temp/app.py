from pages import home, gene, allele, custom_gene, custom_allele, help, citations

data = fetch_data()

data = select_datasets(data) 

if page == pages[0]:
    home.render(data)
elif page == pages[1]:
    gene.render(data)
elif page == pages[2]:
    allele.render(data)
elif page == pages[3]:
    custom_gene.render(data)
elif page == pages[4]:
    custom_allele.render(data)
elif page == pages[5]:
    help.render(data)
elif page == pages[6]:
    citations.render(data)
