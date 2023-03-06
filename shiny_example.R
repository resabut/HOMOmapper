library(shiny)
ui <- fluidPage(
  titlePanel("Title Panel"),
  sidebarLayout(
    sidebarPanel(
      helpText("Sidebar Panel")
    ),
    mainPanel(tabsetPanel(
      tabPanel("tab1",
               fluidRow(
                 column(6,helpText("Col1")),
                 column(6,
                        helpText("Col2"),
                        fluidRow(
                          column(4,style="background-color:#b0c6fb",
                                 helpText("Col1")
                          ),
                          column(4,style="background-color:#ffa153",
                                 helpText("Col2")
                          ),
                          column(4,style="background-color:#b1f6c6",
                                 helpText("Col3")
                          )
                        )
                 )
               )
      ),
      tabPanel("tab2",
               inputPanel(helpText("Input Panel"))
      ),
      tabPanel("tab3",
               wellPanel(helpText("Well Panel"))
      )
    )
    )
  )
)
server <- function(input,output){}
shinyApp(ui=ui,server=server)