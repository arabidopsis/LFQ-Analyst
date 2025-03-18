#!/bin/bash
exec micromamba run -n rlang Rscript \
    -e 'options(shiny.autoload.r = FALSE, shiny.autoreload = TRUE)' \
    -e 'shiny::runApp("." , port=4448, display.mode="normal")'
