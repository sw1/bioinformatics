library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Power"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("alpha",
                  "alpha:",
                  min = .01,
                  max = .25,
                  value = .1)
      ,
      sliderInput("mu",
                  "mean difference:",
                  min = 0,
                  max = 75,
                  value = 25)
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    set.seed(53)
    
    bl <- 200
    md <- input$mu
    
    P1_mu <- bl-md/2
    P2_mu <- bl+md/2
    
    P1 <- rnorm(5000,P1_mu,10)
    P2 <- rnorm(5000,P2_mu,10)
    df <- data.frame(value=c(P1,P2),
                     group=as.factor(rep(c('control','treatment'),each=length(P1))))
    D1 <- density(df %>% filter(group=='control') %>% select(value) %>% unlist())
    D1 <- cbind(D1$x,D1$y)
    D2 <- density(df %>% filter(group=='treatment') %>% select(value) %>% unlist())
    D2 <- cbind(x=D2$x,y=D2$y)
    D <- data.frame(rbind(D1,D2),group=rep(c('control','treatment'),each=nrow(D1)))
    
    alpha <- input$alpha
    
    DF1 <- D %>% filter(group=='control', x > qnorm(1-alpha,P1_mu,10) & x < qnorm(.99999,P1_mu,10))
    DF2 <- D %>% filter(group=='treatment', x > qnorm(0,P1_mu,10) & x < qnorm(1-alpha,P1_mu,10))
    
    p1 <- ggplot(df,aes(value,fill=group)) + 
      geom_density(color='black',alpha=.3) +
      geom_vline(xintercept=qnorm(1-alpha,P1_mu,10),linetype=3) +
      labs(x='',y='',fill='')
    
    if (nrow(DF1)>0) p2 <-  geom_ribbon(data=DF1, aes(x=x,ymin=0,ymax=y),fill='red',alpha=.7) else p2 <- NULL
    if (nrow(DF2)>0) p3 <- geom_ribbon(data=DF2, aes(x=x,ymin=0,ymax=y),fill='blue',alpha=.7) else p3 <- NULL
    
    beta <- sum(D %>% filter(group=='treatment',x > qnorm(0,P2_mu,10) & x < qnorm(1-alpha,P1_mu,10)) %>% select(y))/sum(D %>% filter(group=='treatment') %>% select(y))
    onebeta <- 1-beta
    p4 <- annotate('label',x=235,y=.03,label=paste0('beta==',round(beta,3)),parse=TRUE)
    p5 <- annotate('label',x=235,y=.02,label=paste0('1-beta==',round(onebeta,3)),parse=TRUE)
    
    p1 + p2 + p3 + p4 + p5
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)