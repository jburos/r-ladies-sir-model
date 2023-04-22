#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rstan)
library(tidyverse)
library(geomtextpath)
expose_stan_functions(here::here('model', 'sir-simulation.stan'))

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Simple SIR model simulation"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          sliderInput("N", "Population size:", 
                      min = 100, max = 5000, value = 1000),
          sliderInput("I_0", "Number of initial infected:", 
                      min = 1, max = 50, value = 1),
          sliderInput("beta", "Rate of infection:", 
                      min = 1, max = 10, value = 2, step = 0.1),
          sliderInput("inv_gamma", "Average Recovery Time (days)", 
                      min = 1, max = 50, value = 5, step = 0.1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("sirPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$sirPlot <- renderPlot({
    t <- seq(from = 0.1, to = 50, by = 0.1)
    df <- simulate_sir(t = t, N = input$N, i0 = input$I_0, 
                       beta = input$beta, 
                       gamma = 1/input$inv_gamma) %>%
      map(set_names, c('S', 'I', 'R')) %>% 
      map_dfr(as_tibble_row, .id = 't_index') %>%
      mutate(t_index = as.integer(t_index)) %>%
      left_join(tibble(t = t) %>% mutate(t_index = row_number()))
    df %>%
      tidyr::gather(state, N, S, I, R) %>%
      ggplot(aes(x = t, y = N, colour = state)) +
      geom_labelline(size = 5, aes(label = state)) +
      ggtitle('Simulated Population Dynamics') +
      theme(text = element_text(size = 15))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
