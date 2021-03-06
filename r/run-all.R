# Should all models be run? This will take several minutes
# If FALSE, pre-computed output will be used (recommended)
run <- FALSE

source("r/01_analyze-stomata.R")
source("r/02_analyze-model1.R")
source("r/03_analyze-model2.R")
source("r/04_analyze-model3.R")
source("r/05_plot-deltags.R")

library(rmarkdown)
render("ms/figure-captions.Rmd", word_document(reference_docx = "template.docx"))
render("ms/limonium-si.Rmd", word_document(reference_docx = "template.docx"))
