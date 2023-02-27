# 04_meta-analysis

Meta-analysis pipeline for PACE eczema EWAS. The main meta-analysis report generated can be found here: [meta-analysis-report.html](report/meta-analysis-report.html). The comparison of results across cohorts can be found here: [cohort-comparison-report.html](report/cohort-comparison-report.html).

## Before running analyses

Make sure QC has been run. There should be a QC report for each eczema phenotype:

* [Childhood eczema QC report](../03_qc/report/qc-report-childhood.html)
* [Early-onset eczema QC report](../03_qc/report/qc-report-early_onset.html)
* [Persistent eczema QC report](../03_qc/report/qc-report-persistent.html)

For this analysis, we have chosen to remove the `NEST_W` cohort from model C of the childhood eczema meta-analyses due to inflated beta values and genomic inflation (lambda values). 

To remove cohorts from meta-analyses, remove the cohort results file from the `data/meta_analysis` folder in the RDSF project space.

## Workflow

This uses snakemake on bluecrystal3. Snakemake workflow is presented below:

0. Start a tmux session - `tmux new -s pace-meta`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

`
module add apps/pandoc-2.8.1
snakemake -pr \
-j 9 \
--cluster "qsub \
	-N meta-ewas \
	-l nodes=2:ppn=8,mem=32G,walltime=6:00:00 \
	-o job-errors-and-outputs/meta-ewas-{rule}-error-{jobid} \
	-e job-errors-and-outputs/meta-ewas-{rule}-output-{jobid}" \
--latency-wait 300
`

5. Deactivate the tmux session (CTRL+b + d)
6. Check the completion by going back to the tmux session - `tmux a -t pace-meta`