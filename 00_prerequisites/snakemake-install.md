# Installing Snakemake

This is specific to the high performance computers at University of Bristol. To install it fully elsewhere you'll need to first install [Python](https://www.python.org/) and then [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can follow the instructions, swapping out parts where needs be. 

To check if snakemake is already installed run
`
module add lang/python/anaconda/3.8.5-2021-Jupyter
conda env list
`
and if a snakemake environment (ie. "path_to_env/snakemake") is present then it is installed otherwise continue with the instructions below.

To initiate conda on startup (required for activating environments as specified below) then run `module add lang/python/anaconda/3.8.5-2021-Jupyter; conda init bash` (assuming bash is shell of choice), exit bluepebble and then re-enter.

To insall snakemake, run the following code (change out "/home/tb13101" for your own directory):
`
qsub -I -l select=1:ncpus=1:mem=5G,walltime=1:00:00
module add lang/python/anaconda/3.8.5-2021-Jupyter
cd
conda create -p /user/home/tb13101/conda-envs
conda install -p /user/home/tb13101/conda-envs -c conda-forge mamba
conda activate /user/home/tb13101/conda-envs
mamba create -c conda-forge -c bioconda -n snakemake snakemake
`

To check it's worked run:
`
conda activate /user/home/tb13101/conda-envs/envs/snakemake
snakemake --help
` 
