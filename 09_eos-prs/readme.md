# readme

Before following the snakemake workflow below, you'll need to extract the SNPs and their effects to generate the PRS. The latest PRS of eosinophil counts was made by Xu Y. et al. PMID: [35072137](https://pubmed.ncbi.nlm.nih.gov/35072137/). The results of this PRS (SNPs and weights) can be found in the [PGS catalog](https://www.pgscatalog.org/publication/PGP000051/). These are extracted as part of the workflow (see [`extract-snps.R`](scripts/extract-snps.R)).

snakemake workflow:

0. Start a tmux session - `tmux new -s pace-aries-prs`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:
``` bash
module load apps/plink/1.90
OUTPUT= ## FILL IN
ERROR= ## FILL IN 
snakemake -p \
-j 1 \
--cluster "sbatch \
        --job-name=eos-prs \
        --nodes=1 \
        --ntasks-per-node=1 \
        --cpus-per-task=1 \
        --time=1:00:00 \
        --mem=16GB \
        --output=${OUTPUT}-{rule}-{jobid}.out \
        --error=${ERROR}-{rule}-{jobid}.error"
```
5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t pace-aries-prs`

