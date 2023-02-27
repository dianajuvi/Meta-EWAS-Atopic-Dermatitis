# 03_qc

QC pipeline for PACE eczema EWAS meta-analysis. Reports are generated for each phenotype separately:

* [qc-report-childhood.html](report/qc-report-childhood.html)
* [qc-report-early_onset.html](report/qc-report-early_onset.html)
* [qc-report-persistent.html](report/qc-report-persistent.html)

## Directories setup

Three directories are setup for this analysis: one locally for code editting and some QC, one on RDSF for data storage, one on bluecrystal for the main analyses.

Sort the directory structure in the RDSF PACE data directory so that each cohort has a folder and within these folders are "model" folders with the results in - i.e. m1, m2, m3. Structure example below

```
cohortname
│   README.md    
│	other_files
│
└───m1
│   │   ewas-results-m1a-cohortname.txt.gz
│   │   ewas-results-m1b-cohortname.txt.gz
│   │	ewas-results-m1c-cohortname.txt.gz
│   
└───m2
│   │   ewas-results-m2a-cohortname.txt.gz
│   │   ewas-results-m2b-cohortname.txt.gz
│   │	ewas-results-m2c-cohortname.txt.gz
│
└───m3
    │   ewas-results-m3a-cohortname.txt.gz
    │   ewas-results-m3b-cohortname.txt.gz
    │	ewas-results-m3c-cohortname.txt.gz
    
```

## Workflow

Firstly mannually run through [`manual-data-qc.R`](scripts/manual-data-qc.R) locally after checking through the directory structure on the RDSF. 

Move `conv_file.csv`, a file generated in [`manual-data-qc.R`](scripts/manual-data-qc.R), from the local directory over to the bluecrystal directory

Then the rest of the QC uses snakemake on bluecrystal3. Snakemake workflow is presented below:

0. Start a tmux session - `tmux new -s pace-qc`
1. Activate conda env - `conda activate /newhome/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

`
module add apps/pandoc-2.8.1
snakemake -pr \
-j 9 \
--cluster "qsub \
	-N qc-ewas \
	-l nodes=1:ppn=4,mem=32G,walltime=4:00:00 \
	-o job-errors-and-outputs/qc-ewas-{rule}-error-{jobid} \
	-e job-errors-and-outputs/qc-ewas-{rule}-output-{jobid}"
`

5. Detach the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t pace-qc`
