library(shiny)
library(toxEval)
library(shinydashboard)
library(shinyWidgets)
library(tidyverse)
library(DT)

header <- dashboardHeader(
    title = "Explore Mixtures"
)

sidebar <- dashboardSidebar(
    sidebarMenu(
        numericInput("n_sites",
                     "Site Threshold:",
                     min = 1,
                     max = 75,
                     value = 10),
        numericInput("ear_thresh",
                      label = "EAR threshold:",
                      min = 0.000001,
                      max = 0.001,
                     value = 0.0001)
    )
)

body <- dashboardBody(
    tabBox(width = 12, id="mainOut",
           tabPanel(title = "ToxCast", value = "toxCast",
                    DT::dataTableOutput(outputId = "epMixes"),
                    plotOutput("tox_missing",
                               width = "100%",
                               height = "750px")
           ),
           tabPanel(title = "AOP", value = "AOP",
                    DT::dataTableOutput(outputId = "aopMixes"),
                    plotOutput("aop_missing",
                               width = "100%",
                               height = "300px")                    
           ),
           tabPanel(title = "GeneID", value = "gene",
                    DT::dataTableOutput(outputId = "geneMixes"),
                    plotOutput("gene_missing",
                               width = "100%",
                               height = "450px")                      
           ),
           tabPanel(title = "Panther", value = "panther",
                    DT::dataTableOutput(outputId = "pantherMixes"),
                    plotOutput("panther_missing",
                               width = "100%",
                               height = "450px")                      
           )
    )
)

dashboardPage(header, sidebar, body)
