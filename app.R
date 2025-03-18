VERSION <- "v1.5"

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

`%>%` <- magrittr::`%>%`

source("R/functions.R")
source("R/volcano_function.R")
source("R/tests.R")
source("R/enrichment_functions.R")
source("R/lfq_ui.R")
source("R/lfq_server.R")




server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024^2) ## Set maximum upload size to 100MB
  moduleServer("lfq", lfq_server)
  # server_bg(input, output, session)
}

app <- shinyApp(ui = lfq_ui(), server = server)
