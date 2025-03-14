VERSION <- "v1.3"

library("shiny")
loadNamespace("DEP")
# library("SummarizedExperiment")
# library("tidyverse")
# library("DEP")
# library("testthat")
# library("shinydashboard")
# library("shinyjs")
# library("shinyalert")
# library("ComplexHeatmap")
# library("limma")
# library("DT")
# library("ggrepel")
# library("httr")
# library("rjson")
# library("svglite")
# library("shinycssloaders")
# library("shiny.info")

source("R/functions.R")
source("R/volcano_function.R")
source("R/tests.R")
source("R/demo_functions.R")
source("R/enrichment_functions.R")


`%>%` <- magrittr::`%>%`