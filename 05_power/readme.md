# 06_power

The scripts and reports in this directory aim to assess the power to detect associations between DNAm and the three AD subtypes using the data we have.

The supplementary note for the power calculations can be found in [supplementary-note.pdf](report/supplementary-note.pdf).

To generate the supplementary note run the following:

```bash
cd report
Rscript -e "rmarkdown::render('supplementary-note.Rmd', params = list(study_data = 'report_data/study_data.tsv', power_levels = 'report_data/power-levels.tsv'))"
``