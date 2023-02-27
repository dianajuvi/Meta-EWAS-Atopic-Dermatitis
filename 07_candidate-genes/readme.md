# Candidate gene analysis readme

In this repo we assess whether DNA methylation surrounding candidate eczema genes identified by [Sobczyk M et al.](https://www.sciencedirect.com/science/article/pii/S0022202X2101160X) associate with eczema more strongly than DNA methylation elsewhere in the genome.

## Before running analyses

Make sure the meta-analysis has been run. There should be a [report](../04_meta-analysis/report/meta-analysis-report.html) and there should be results in `04_meta-analysis/results/metal-res/` for each eczema phenotype and model.

## Workflow

This uses snakemake on bluepebble. Snakemake workflow is presented below:

0. Start a tmux session - `tmux new -s pace-cg`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

``` bash
module add apps/pandoc-2.8.1
OUTPUT= ## FILL IN
ERROR= ## FILL IN 
snakemake -pr \
-j 1 \
--cluster "sbatch \
#	--test-only \
	--job-name=candidate-genes \
	--nodes=1 \
	--ntasks-per-node=1 \
	--cpus-per-task=1 \
	--time=1:00:00 \
	--mem=16GB \
	--output=${OUTPUT}-{rule}-{jobid}.out \
	--error=${ERROR}-{rule}-{jobid}.error" \
--latency-wait 300
```

5. Deactivate the tmux session (CTRL+b + d)
6. Check the completion by going back to the tmux session - `tmux a -t pace-cg`
