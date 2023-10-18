
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(data.table)
library(ComplexUpset)

ui <- dashboardPage(
  dashboardHeader(title = "FOXP2 Analysis"),
  
  dashboardSidebar(  # Sidebar definition
    # sidebarMenu(
    #   menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    #   menuItem("Settings", icon = icon("gear")),
    #   menuItem("Help", icon = icon("question-circle"))
    # ),
    radioButtons("selectOpt", "Select dataset:",
                 list("Evolution"='evolution', "Time"='time'),
                 selected = 'evolution'),
    
    conditionalPanel(
      condition = "input.selectOpt=='evolution'",
      selectizeInput(
        inputId = "data_evo", 
        label = "Select data evolution",
        choices = c("Chimp_vs_Control","Human_vs_Control","KE_vs_Control","Mouse_vs_Control",   
                    "N303T_vs_Control","PhyDis_vs_Control","RouAeg_vs_Control","S325N_vs_Control",   
                    "Songbird_vs_Control"),
        multiple = T,
        selected = c("Chimp_vs_Control","Human_vs_Control","KE_vs_Control","Mouse_vs_Control")
      )
    ),
    conditionalPanel(
      condition = "input.selectOpt=='time'",
      selectizeInput(
        inputId = "data_time", 
        label = "Select data time",
        choices = c("foxp2-diff_vs_undif","ke-diff_vs_undif","ctrl-diff_vs_undif","foxp2_ctr_undif",    
                    "ke_ctr_undif", "foxp2_ctr_diff","foxp2_ke_undif","ke_ctr_diff",        
                    "foxp2_ke_diff"),
        multiple = T,
        selected = c("foxp2-diff_vs_undif","ke-diff_vs_undif")
      )
    ),
    radioButtons(
      inputId = "booleanCounts", 
      label = ("Display counts:"), 
      choices = c("TRUE","FALSE"),
      selected = "FALSE"
    ),
    selectInput(
      inputId = "sortBy", 
      label = ("Sort by:"), 
      choices = c("degree", "cardinality"), 
      multiple = F,
      selected = "degree"
    ),
    sliderInput(
      inputId = "nSets",
      label = ("Min. sets:"),
      min = 50,
      max = 500,
      value = 50,
      step = 50,
      round = F,
      ticks = T,
      animate = T,
      width = NULL
    ),
    sliderInput(
      inputId = "angle", 
      label = ("Font angle:"),
      min = 0,
      max = 90,
      value = 0,
      step = 15)
  ),
  
  dashboardBody(
    # Top box spanning entire width
    fluidRow(
      box(title = "Upset Plot", width = 12, status = "primary", "", plotOutput(outputId = "Plot") ),
      tags$p(
        downloadButton(
          outputId = "downloadPlot",
          label = "Download"),
        align="right")
    ),
    # Two boxes below
    fluidRow(
      box(title = "Table with set intersections", width = 12, status = "warning", "", DTOutput("table") ),
      hr(),
      tags$p(
        downloadButton(
          outputId = "downloadTable",
          label = "Download"),
        align="right")
    )
  )
)

# Server
server <- function(input, output) {
  
  f_read_data <- function(dat) {
    d <- readRDS( file.path( "data", paste0("dge_", dat, ".rds") ) )
    dt <- as.data.table( d, key = "gene_id" )
    dt_w <- dcast( dt, formula = gene_name + direction ~ group, fun.aggregate = length, value.var = "direction" )
    return(list(dt, dt_w))
  }  
  
  setsSelected <- reactive({
    if(input$selectOpt == "evolution") {
      input$data_evo
    } else {
      input$data_time
    }
  })
  
  dataSelected <- reactive({
    if(input$selectOpt == "evolution") {
      f_read_data("evolution")[[1]]
    } else {
      f_read_data("time")[[1]]
    }
  })
  
  dataSelected_wide <- reactive({
    if(input$selectOpt == "evolution") {
      f_read_data("evolution")[[2]]
    } else {
      f_read_data("time")[[2]]
    }
  })
  
  # upSetPlot ------------------------------------------------------------------ #
  
  plotInput <- shiny::reactive({ ComplexUpset::upset(
    data = dataSelected_wide(),
    intersect = setsSelected(),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts=input$booleanCounts, mode = "intersect", mapping=aes(fill=direction),
        text=list(
          angle=input$angle, 
          size=3
        )
      )
      + theme(
        axis.ticks.y=element_line(),
        axis.title.y=element_text(face = "bold"),
        legend.title=element_text(face = "bold")
      )
      + ylab('# deregulated genes')
      + guides(fill=guide_legend(title='direction'))
      + scale_fill_manual(values = c("DOWN" = "darkgreen", "UP" = "red")) 
      + geom_hline(yintercept = seq(50, 1250, by=250),
                   linetype = 2, colour = "gray", show.legend = TRUE)
    ),
    name = "",
    annotations = list(),
    themes = upset_modify_themes(
      list(
        'intersections_matrix'=theme(text=element_text(size=15)
        ))
    ),
    labeller = identity,
    height_ratio = 0.5,
    width_ratio = 0.25,
    wrap = FALSE,
    set_sizes=(
      upset_set_size()
      + theme(
        axis.ticks.x=element_line(),
        axis.title.x=element_text(face = "bold"),
        text=element_text(face = "bold"))
      + ylab('set size')
    ),
    mode = "intersect",
    queries = list( upset_query(intersect=setsSelected(), color='orange') ),
    guides = NULL,
    encode_sets = TRUE,
    min_size = input$nSets,
    sort_intersections_by=input$sortBy,
    sort_sets='descending'
  )
  })
  
  output$Plot <- shiny::renderPlot({
    print(plotInput())
  })
  
  # DT table reactive ---------------------------------------------------------- #
  f_dt_tbl <- function(dat_wid, dat_ori,sset) {
    s <- sset
    s_len <- length( s )
    s_query <- dat_wid[rowSums(dat_wid[, ..s]) == s_len][,1]
    out <- dat_ori[dat_ori[, gene_name %in% s_query[[1]] & group %in% s]]
    setDF( out )
    return(out)
  }
  
  tableInput <- shiny::reactive({
    res = f_dt_tbl( dataSelected_wide(), dataSelected(), setsSelected() )
  })
  
  output$table <- renderDT({
    datatable( tableInput() ) %>% formatRound(4:7, 3) %>% formatSignif(8:9, 3) %>%
      formatStyle(
        'direction',
        backgroundColor = styleEqual(c("DOWN", "UP"), c('green', 'red')) )
  }, options = list(
    pageLength = 10,
    lengthMenu = c("10","25","50","100"),
    initComplete = I("function(settings, json) {alert('Done.');}")
  )
  )
  
  # Download Plot -------------------------------------------------------------- #
  output$downloadPlot <- shiny::downloadHandler(
    filename = function(){
      paste0("UpSetPlot-", Sys.time(), ".pdf")
    },
    contentType = "image/pdf",
    content = function(file){
      pdf(file, onefile = F, width = 16, height = 8)
      print(plotInput())
      dev.off()
    })
  
  # Download Table ------------------------------------------------------------- #
  output$downloadTable <- shiny::downloadHandler(
    filename = function(){
      paste0("UpSetDataSet-", Sys.time(), ".csv")
    },
    content = function(file){
      write.csv(tableInput(), file, row.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)