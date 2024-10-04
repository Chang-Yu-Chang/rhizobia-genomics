import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

# Config file and sample sheets
configfile: "../config/config.yaml"
#validate(config, schema = "../schemas/config.schema.yaml")

genome_ids = pd.read_csv(config["samples"])["genome_id"].tolist()
#print(genome_ids)
#validate(samples, schema = "../schemas/samples.schema.yaml")
