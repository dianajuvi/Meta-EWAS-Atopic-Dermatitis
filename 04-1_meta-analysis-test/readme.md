# 04-1_meta-analysis-test

Meta-analysis pipeline for PACE eczema EWAS. Some cohorts defined eczema using doctor diagnosis, others used rash symptoms and others used a combination of both. The code in this directory is designed to examine the differences in EWAS meta-analysis for cohorts who used doctor-only diagnoses and the other cohorts.

The report generated can be found here: [meta-analysis-report.html](report/meta-analysis-comparison-report.html).

`
snakemake -prk \
-j 2 \
--cluster "sbatch \
  --job-name=pace-meta \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task=1 \
  --time=2:00:00 \
  --mem=32G" 
`

Steps:
1. UPDATE THE COHORTS IN EACH OF THE RDSF FOLDERS TO MAKE SURE WE HAVE A DOCTOR-DIAGNOSIS ONLY META!!
2. Re-run it 
3. Check report
4. Add in persistent and early-onset