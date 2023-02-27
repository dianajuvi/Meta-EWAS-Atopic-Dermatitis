# 04-2_meta-analysis-test2

Meta-analysis pipeline for PACE eczema EWAS. Some cohorts defined eczema using the Hanifin-Rajka criteria, which is thought by some to be a superior definition of eczema case status compared to other definitions (including doctor diagnosis). The cohorts that used this criteria are: EDEN, GOYA, IOW-F1, IOW-F2. The code in this directory is designed to examine the differences in EWAS meta-analysis for cohorts who used the Hanifin-Rajka criteria and the other cohorts.

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