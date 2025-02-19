conda activate tobias_snakemake_env
snakemake --configfile config.yaml --use-conda --cores 30 --conda-prefix Temp --keep-going
