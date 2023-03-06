library(shiny)
library(shinyBS)
ui <- fluidPage(titlePanel("ROH calculator"),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Input files",
                                           helpText("Please select the input file format."),
                                           tabsetPanel(
                                             tabPanel(
                                               "Binary format (.bed)",
                                               
                                               fileInput("bed_file", "Choose .bed file",
                                                         accept = ".bed"),
                                               fileInput("bim_file",
                                                         "Choose .bim file",
                                                         accept = ".bim"),
                                               fileInput("fam_file",
                                                         "Choose .fam file",
                                                         accept = ".fam")
                                             ),
                                             tabPanel(
                                               "PLINK format (.ped)",
                                               fileInput("ped_file", "Choose .ped file",
                                                         accept = ".ped"),
                                               fileInput("map_file", "Choose .map file",
                                                         accept = ".map")
                                             )
                                           )),
                                  tabPanel("Configuration",
                                           inputPanel(helpText("Input Panel")),
                                           
                                           bsButton(
                                             inputId = "StartRunButtonlist(",
                                             label = "Run",
                                             icon = icon("play"), # Optional
                                             style = "default",
                                             size = "large",
                                             disabled = TRUE
                                           )),
                                  tabPanel("tab3",
                                           wellPanel(helpText("Well Panel")))
                                )
                              ))
server <- function(input, output) {
}
shinyApp(ui = ui, server = server) 