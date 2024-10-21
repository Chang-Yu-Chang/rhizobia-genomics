
from unipressed import IdMappingClient

request = IdMappingClient.submit(
    source = "Gene_Name", 
    dest = "UniProtKB", 
    ids = {"COQ3"},
    taxon_id = "266834"
)
request
list(request.each_result())


import pandas as pd

# Example list of data as a list of dictionaries
data = [
    {"Gene": "TP53", "UniProt ID": "P04637"},
    {"Gene": "BRCA1", "UniProt ID": "P38398"},
    {"Gene": "EGFR", "UniProt ID": "P00533"}
]

# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame(data)

# Specify the file name and write to CSV
df.to_csv('genes.csv', index=False)

print('Data written to genes.csv')
