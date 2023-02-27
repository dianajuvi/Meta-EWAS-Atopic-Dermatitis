# readme

Before following the snakemake workflow below, you'll need to extract the SNPs and their effects to generate the PRS. These can be extracted using [`extract-snps.R`](scripts/extract-snps.R)

snakemake workflow:

0. Start a tmux session - `tmux new -s pace-aries-prs`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:
`
module add apps/pandoc-2.8.1
module load apps/plink-1.90
module add apps/plink-2.00
snakemake -p \
-j 3 \
--cluster "qsub \
	-N pace-aries-prs \
	-l nodes=1:ppn=6,mem=32G,walltime=8:00:00 \
	-o job-errors-and-outputs/pace-aries-prs-{rule}-error-{jobid} \
	-e job-errors-and-outputs/pace-aries-prs-{rule}-output-{jobid}"
`
5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t pace-aries-prs`

