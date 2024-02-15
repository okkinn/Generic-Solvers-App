# John Nico T. De Castro
# 2022-12523
# B1L     Project
# ui side

library(shiny)
library(shinydashboard)
library(DT)
library(shinyWidgets)

foodChoices = read.csv('nutritional_values.csv')["Foods"]

dashboardPage(
  dashboardHeader(title = "Generic Solvers"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName="A"),
      menuItem("Quadratic Spline Interpolation", tabName="QSI"),
      menuItem("Polynomial Regression", tabName = "PR"),
      menuItem("Diet Problem Solver", tabName = "DPS")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="A",
        h2("About", style="text-align:center"),
        br(),
        div("The application is designed to provide users with generic solvers for quadratic spline interpolation and polynomial regression, and a simplex implementation for the diet optimization problem. It consists of an RShiny application and a shinyApp web application.", style="text-align:center"),
        br(),
        div("Developed by John Nico T. De Castro as a project requirement for CMSC 150.", style="text-align:center"),
        br(),
        div("You may access the user manual at this link: https://tinyurl.com/DeCastro150ProjectUserManual", style="text-align:center")
      ),
      
      tabItem(tabName="QSI",
        sidebarLayout(
          sidebarPanel(
            fileInput("QSIfile", "Choose CSV File", accept = c(".csv")),
            numericInput("QSIx_estimate", "X Value to Estimate", value = 0),
            sliderInput("QSIround", "Round Value by Corresponding Decimals", min = 1, max = 10, value = 6),
            actionButton("QSIcalculate", "Calculate"),
            width = 3
          ),
          
          mainPanel(
            box(dataTableOutput("QSItable")),
            box(htmlOutput("QSIintro_text"), htmlOutput("QSIintervals_output"), textOutput("QSIestimate_output"), textOutput("QSIfunction_output_text"), textOutput("QSIfunction_output")),
            box(plotOutput("QSIplot"))
          )
        )
      ),
      
      tabItem(tabName="PR",
        sidebarLayout(
          sidebarPanel(
            fileInput("PRfile", "Choose CSV File", accept = c(".csv")),
            sliderInput("PRdegree", "Polynomial Degree", min = 1, max = 10, value = 4),
            numericInput("PRx_estimate", "X Value to Estimate", value = 0),
            sliderInput("PRround", "Round Value by Corresponding Decimals", min = 1, max = 10, value = 2),
            actionButton("PRcalculate", "Calculate"),
            width = 3
          ),
          
          mainPanel(
            box(dataTableOutput("PRtable")),
            box(htmlOutput("PRfunction_output"), textOutput("PRestimate_output")),
            box(plotOutput("PRplot"))
          )
        )
      ),
       
      tabItem(tabName="DPS",
        fluidRow(
          column(width=6,
            box(
              pickerInput(
                inputId = "foodInput", label = "Select Foods", choices = read.csv('nutritional_values.csv')["Foods"],
                selected = read.csv('nutritional_values.csv')[c(1,2,3,4,5,9,7,8,6,10,11,12,13,14,15,16,17,18,19,20),1],
                multiple = TRUE,
                options = pickerOptions(actionsBox = TRUE, size=10, selectedTextFormat = "count > 3")
              ),
              actionButton("DPScalculate", "Start Solve"),
              width=12
            ),
            box(dataTableOutput("FoodTable"), width=12),
          ),
          column(width=6,
                 box(textOutput("DPSzVal"), width=12),
                 box(dataTableOutput("DPStable"), width=12)      
          )
        ),
        fluidRow(
          uiOutput("panels") # initialtableau, tableaus, basic solutions, final solution
        )
      )
    )
  )
)