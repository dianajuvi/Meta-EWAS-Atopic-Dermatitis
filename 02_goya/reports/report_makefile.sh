#!/bin/bash

# -------------------------------------------------------
# Report makefile
# -------------------------------------------------------

cd ~/PACE/pace_ad/reports/

RMDFILE=pace_ad_report

Rscript -e "require(knitr); require(markdown); knit('${RMDFILE}.Rmd', '${RMDFILE}.md'); markdownToHTML('${RMDFILE}.md', '${RMDFILE}.html', options=c('use_xhtml', 'base64_images'))"
