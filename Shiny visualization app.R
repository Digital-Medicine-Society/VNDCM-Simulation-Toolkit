#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)

datasets <- c("economics", "faithfuld", "seals")
# Define UI for application
# ui <- page_sidebar(
#   sidebar = sidebar(
#     varSelectInput("xvar", "X variable", df_num),
#     varSelectInput("yvar", "Y variable", df_num),
#     plotOutput("scatter")
#   )
#     
#   server <- function(input, output, session) {
#     dataset <- reactive({
#           get(input$dataset, "package:ggplot2")
#         })
#   }
  
ui <- fluidPage(
  selectInput("dataset", "Dataset", choices = datasets),
  verbatimTextOutput("summary"),
  plotOutput("plot")
)

server <- function(input, output, session) {
  dataset <- reactive({
    get(input$dataset, "package:ggplot2")
  })
  output$summary <- renderPrint({
    summary(dataset())
  })
  output$plot <- renderPlot({
    ggplot(dataset())
  }, res = 96)
}

shinyApp(ui, server)


# ui <- fluidPage(
#   sliderInput("x", label = "If x is", min = 1, max = 50, value = 30),
#   "then x times 5 is",
#   textOutput("product")
# )
# 
# server <- function(input, output, session) {
#   output$product <- renderText({
#     input$x * 5
#   })
# }
# 
# shinyApp(ui, server)

# ui <- fluidPage(
#  
#   textInput("name", "What's your name?"),
#   textOutput("greeting")
#   
# )
# 
# server <- function(input, output, session) {
#   
#   output$greeting <- renderText({
#     paste0("Hello ", input$name)
#   })
#   
# }
# 
# 
# shinyApp(ui, server)
