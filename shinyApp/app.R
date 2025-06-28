library(shiny)
library(poem)
library(ggplot2)
library(shinycssloaders)
library(shinydashboard)
library(shinyjqui)
library(BiocParallel)
library(DT)


################################## INITIATION

# base plot fontsize
base_size <- 15

# colors
cols <- c("#F21A00"=1L, "#EBCC2A"=2L, "#78B7C5"=3L)
colChoices <- c("Red"=1L, "Yellow"=2L, "Blue"=3L)

# starting matrices (should be squared and all of the same dims!)
m2 <- m3 <- m <- matrix(1L, nrow=10, ncol=10)
m3[lower.tri(m)] <- m[lower.tri(m)] <- 2L
m2[lower.tri(m,diag = TRUE)] <- 2L
m3[3:4,1:7] <- 2L
m3[4,8] <- 2L

data(metric_info)
metric_info <- metric_info[,c(1,4)]


################################## REQUIRED FUNCTIONS

getStats <- function(mats){
  x <- lapply(mats, \(x) reshape2::melt(x, varnames=c("y","x")))
  x <- as.data.frame(cbind(x[[1]][,1:2],sapply(x, \(x) as.factor(x[,3]))))
  res1 <- dplyr::bind_rows(list(
    p1=poem::getPartitionMetrics(x$h1, x$h2, level = "class"),
    p2=poem::getPartitionMetrics(x$h1, x$h3, level = "class")), 
    .id="Partition")
  res1$cl <- apply(res1[,c("class","cluster")], 1, \(x) x[!is.na(x)])
  res1$class <- res1$cluster <- NULL
  res1 <- reshape2::melt(res1, id.vars=c("Partition","cl"))
  res1b <- dplyr::bind_rows(list(
    p1=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h2, location=x[,2:1],
                                       level = "class", verbose=FALSE),
    p2=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h3, location=x[,2:1],
                                       level = "class", verbose=FALSE)), 
    .id="Partition")
  res1b$cl <- apply(res1b[,c("class","cluster")], 1, \(x) x[!is.na(x)])
  res1b$class <- res1b$cluster <- NULL
  res1b <- reshape2::melt(res1b, id.vars=c("Partition","cl"))
  res1 <- rbind(res1,res1b)
  res1 <- res1[!is.na(res1$value),]
  
  res2 <- dplyr::bind_rows(list(
    p1=poem::getPartitionMetrics(x$h1, x$h2, level = "dataset", metrics=c("RI","ARI")),
    p2=poem::getPartitionMetrics(x$h1, x$h3, level = "dataset", metrics=c("RI","ARI"))),
    .id="Partition")
  res2 <- cbind(res2[,-1], dplyr::bind_rows(list(
    p1=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h2, location=x[,2:1], level = "dataset",
                                       metrics=c("nsRI", "nsARI", "SpatialRI", "SpatialARI"), verbose=FALSE),
    p2=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h3, location=x[,2:1], level = "dataset",
                                       metrics=c("nsRI", "nsARI", "SpatialRI", "SpatialARI"), verbose=FALSE)), 
    .id="Partition"))
  res2$cl <- "global"
  #res2 <- reshape2::melt(res2, id.vars=c("Partition","cl"))
  
  res3 <- dplyr::bind_rows(list(
    p1=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h2, location=x[,2:1], level = "dataset", original=TRUE,
                                       metrics=c("SpatialRI", "SpatialARI"), verbose=FALSE),
    p2=poem::getSpatialExternalMetrics(true=x$h1, pred=x$h3, location=x[,2:1], level = "dataset", original=TRUE,
                                       metrics=c("SpatialRI", "SpatialARI"), verbose=FALSE)), 
    .id="Partition")
  colnames(res3)[2:3] <- paste0(colnames(res3)[2:3], "(original)")
  res3$cl <- "global"
  rbind(res1,reshape2::melt(cbind(res2,res3[,2:3]), id.vars=c("Partition","cl")))
}

getSpotwise <- function(mats){
  y <- lapply(mats, \(x) reshape2::melt(t(x), varnames=c("y","x")))
  
  x <- lapply(y[-1], \(x){
    x$SPC <- poem:::getPairConcordance(y[[1]]$value, x$value, useNegatives = TRUE)
    tmp <- getSpatialExternalMetrics(true=y[[1]]$value, pred=x$value,
                                     location=x[,1:2], level="element", 
                                     metrics=c("nsSPC"))
    x$nsSPC <- tmp[,ncol(tmp)]
    x$spatialSPC <- spatialARI(y[[1]]$value, x$value, x[,1:2], spotWise=TRUE)
    x[["spatialSPC(original)"]] <- spatialARI(y[[1]]$value, x$value, x[,1:2],
                                              spotWise=TRUE, original=TRUE)
    x
  })
  names(x) <- c("Partition p1", "Partition p2")
  x <- dplyr::bind_rows(x, .id="Partition")
  reshape2::melt(x[,-4], id.vars=c("Partition","y","x"))
}

################################## UI

ui <- fluidPage(
  titlePanel("poem - metrics explorer"),
  tags$head(tags$style(("
"))),
  mainPanel(style="width: 100%;", 
    fluidRow(
      column(4, style="text-align: center;", 
        tags$h3("Ground truth"),
        plotOutput("heatmap1", click="heatmap_click1", height = "300px",
                   brush=brushOpts("heatmap_brush1", resetOnNew=TRUE, fill=NA))),
      column(4, style="text-align: center;", 
        tags$h3("Partition p1"),
        plotOutput("heatmap2", click="heatmap_click2", height = "300px",
                   brush=brushOpts("heatmap_brush2", resetOnNew=TRUE, fill=NA))),
      column(4, style="text-align: center;", 
        tags$h3("Partition p2"),
        plotOutput("heatmap3", click="heatmap_click3", height = "300px",
                   brush=brushOpts("heatmap_brush3", resetOnNew=TRUE, fill=NA)))
    ),
    fluidRow(style="margin-top: 20px;",
      column(6, selectizeInput("newColor", label="Click (or drag) on images to set to :", choices=colChoices)),
      column(6, style="text-align: right;", actionButton("reset", "Reset"))),
    #tags$hr(style="width: 80%; color: #000000; border: 1px solid black;"),
    tabsetPanel(
      tabPanel("Dataset- and Class-level metrics", 
               fluidRow(column(4, actionButton("compile", "Compute metrics")),
                        column(6, style="text-align: right;", checkboxInput("free","Free y-axis", value=FALSE))),
               withSpinner(jqui_resizable(plotOutput("statsPlot", height="600px")))),
      tabPanel("Element-level metrics", 
               fluidRow(column(4, actionButton("compileSpotwise", "Compute metrics")),
                        column(6, style="text-align: right;", checkboxInput("labs","Write values", value=TRUE))),
               withSpinner(withSpinner(plotOutput("spotwisePlot", height="600px")))),
      tabPanel("Metric descriptions", DTOutput("metricsInfo")),
      tabPanel("About", tags$div(style="padding: 20px;",
        tags$p("This small app is meant to help understand the behavior of the (spatial and non-spatial) partition metrics implemented in the ",
               tags$a(href="https://bioconductor.org/packages/poem/", "poem"), " bioconductor package."),
        tags$p("For more in-depth discussion of the metrics, see the ",
               tags$a(href="https://bioconductor.org/packages/release/bioc/vignettes/poem/inst/doc/poem.html", "package vignettes"),
               " and the ", tags$a(href="https://www.biorxiv.org/content/10.1101/2024.11.28.625845v1", "accompanying manuscript"), ".")
              )),
    )
  )
)

################################## SERVER

server <- function(input, output, session) {
  
  mats <- reactiveValues(h1=m, h2=m2, h3=m3)
  stats <- reactiveVal(NULL)
  spotwise <- reactiveVal(NULL)
  helpShown <- reactiveVal(FALSE)

  hmClickEvent <- function(event, name){
    if(is.null(event)) return(NULL)
    m <- mats[[name]]
    pos <- round(c(event$y, event$x))
    if(m[pos[2], pos[1]] != as.integer(input$newColor)){
      m[pos[2], pos[1]] <- as.integer(input$newColor)
      mats[[name]] <- m
      stats(NULL)
      spotwise(NULL)
    }
  }

  hmBrushEvent <- function(event, name){
    if(is.null(event)) return(NULL)
    stats(NULL)
    spotwise(NULL)
    pos <- round(c(event$xmin, event$xmax, event$ymin, event$ymax))
    m <- mats[[name]]
    m[pos[1]:pos[2], pos[3]:pos[4]] <- as.integer(input$newColor)
    mats[[name]] <- m
  }

  drawHm <- function(name, title=NULL){
    m <- mats[[name]]
    par(mar=rep(0,4))
    graphics::image(x=seq_len(ncol(m)), y=seq_len(nrow(m)), z=m, #main=title, 
          col=names(cols), breaks=c(0,1.5,2.5,3.5), asp=1, ylab=NA, xlab=NA,
          xaxt="n", yaxt="n", bty="n")
  }
  
  output$heatmap1 <- renderPlot({ drawHm("h1") })
  output$heatmap2 <- renderPlot({ drawHm("h2") })
  output$heatmap3 <- renderPlot({ drawHm("h3") })
  
  output$statsPlot <- renderPlot({
    s <- stats()
    if(is.null(s)) return(NULL)
    col2 <- setNames(names(cols), as.character(as.integer(cols)))
    ggplot(s, aes(Partition, value, colour=cl, group=cl)) + geom_point() + 
      geom_line() + facet_wrap(~variable, scales=ifelse(input$free, "free_y", "fixed")) + 
      scale_color_manual(values=c("global"="black", col2)) + 
      scale_y_continuous(breaks=scales::pretty_breaks(3)) +
      theme_bw(base_size = base_size)
  })
  
  output$spotwisePlot <- renderPlot({
    y <- spotwise()
    if(is.null(y)) return(NULL)
    y$value <- 100*y$value
    p <- ggplot(y, aes(x,y,fill=value,label=round(value))) + geom_tile()
    if(input$labs) p <- p + geom_text()
    p + facet_grid(Partition~variable) + theme_void(base_size) + 
      scale_fill_viridis_c(option="C", direction=1, breaks=scales::pretty_breaks(4)) +
      labs(fill="value (%)") + theme(strip.text.y=element_text(angle=90))
  })
  
  output$metricsInfo <- renderDT({
    metric_info
  })

  observeEvent(input$heatmap_click1, hmClickEvent(input$heatmap_click1, "h1"))
  observeEvent(input$heatmap_click2, hmClickEvent(input$heatmap_click2, "h2"))
  observeEvent(input$heatmap_click3, hmClickEvent(input$heatmap_click3, "h3"))
  observeEvent(input$heatmap_brush1, hmBrushEvent(input$heatmap_brush1, "h1"))
  observeEvent(input$heatmap_brush2, hmBrushEvent(input$heatmap_brush2, "h2"))
  observeEvent(input$heatmap_brush3, hmBrushEvent(input$heatmap_brush3, "h3"))
  
  observeEvent(input$reset, {
    mats[["h1"]] <- m
    mats[["h2"]] <- m2
    mats[["h3"]] <- m3
    stats(getStats(reactiveValuesToList(mats)))
    spotwise(getSpotwise(reactiveValuesToList(mats)))
  })
  
  observeEvent(input$compile, {
    stats(getStats(reactiveValuesToList(mats)))
  })
  observeEvent(input$compileSpotwise, {
    spotwise(getSpotwise(reactiveValuesToList(mats)))
  })
}

shinyApp(ui = ui, server = server)
