# 01_aries

EWAS of eczema in ARIES. Workflow seen below. Output can be seen in [ewas-report.html](report/ewas-report.html).

snakemake workflow:

0. Start a tmux session - `tmux new -s aries-ewas`
1. Activate conda env - `conda activate /newhome/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -n`
4. Submit pipeline as a job:
`
module add apps/pandoc-2.8.1
snakemake -p \
-j 9 \
--cluster "qsub \
	-N aries-ewas \
	-l nodes=1:ppn=6,mem=32G,walltime=8:00:00 \
	-o job-errors-and-outputs/aries-ewas-{rule}-error-{jobid} \
	-e job-errors-and-outputs/aries-ewas-{rule}-output-{jobid}"
`
5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t aries-ewas`

Problem: want to add date to EWAS output scripts, but this would cause snakemake to re-run the EWAS each time even though nothing else has changed... Potential solution could involve having another rule to change the name to add the date if the results files have changed.

