# app.R — application entry point (launched via shiny::runApp on this directory)
# In the lab/HPC environment, .libPaths() and setwd() were set here to the project library and app directory, but removed for release. Run from the app dir.

library(shiny)            # interactive web application framework
library(Seurat)           # single-cell / spatial analysis and visualisation
library(ggplot2)          # layered grammar-of-graphics plotting
library(patchwork)        # compose multiple ggplots into one layout
library(plotly)           # interactive web-based graphs
library(DT)               # interactive DataTables in Shiny
library(dplyr)            # data manipulation
library(scCustomize)      # single-cell visualisation utilities
library(shinycssloaders)  # loading spinners for outputs
library(shinyjs)          # JavaScript helpers within Shiny
library(future)           # parallel / asynchronous computation
library(ggnewscale)       # multiple colour/fill scales in one ggplot
library(BPCells)          # on-disk compressed matrices for large datasets
library(harmony)          # batch correction / integration
library(Azimuth)          # reference-based cell type annotation
library(clusterProfiler)  # functional enrichment analysis
library(org.Mm.eg.db)     # mouse gene annotation database

options(shiny.maxRequestSize = 10000 * 1024^2)  # allow large (~10 GB) uploads

# Source shared globals and feature modules (each separately redacted)
source("global.R")
source("queue_module.R")
source("qc_module.R")
source("analysis_module.R")
source("export_module.R")
source("azimuth_module.R")
source("multi_sample_module.R")
source("annotation.R")
source("ui_server.R")

ui     <- app_ui()
server <- app_server
shinyApp(ui, server)
