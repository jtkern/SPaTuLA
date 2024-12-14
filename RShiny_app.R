#capstone app 12/3/24
#inclues full SPaTuLA model

#libraries required
library(shiny)
library(shinyMatrix)
library(colourpicker)
library(shinydashboard)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(readr)
library(vembedr)
library(bslib)
library(tibble)

#model source code
source('~/KaitlynRRStudio/BIOE_Capstone/SPaTuLA_functions.R')

# Define UI for application
ui <- dashboardPage(
  
  skin="black",
  
  # Application title and pages defined
  dashboardHeader(title="SPaTuLA"),
  dashboardSidebar(  
    sidebarMenu(
      id="tabs",
      menuItem("About", tabName="About", icon=icon("user")),
      menuItem("Tutorial", tabName="Tutorial", icon=icon("dashboard")),
      menuItem("SPaTuLA Simulator", tabName="Simulator", icon=icon("table"), startExpanded = TRUE,
               menuSubItem("Inputs", tabName = "Inputs"),
               menuSubItem("Outputs", tabName = "Outputs"))
    )
  ),
  
  #defines contents and aesthetics of each page (aka tab)
  dashboardBody(
    tabItems(
      
      #about page
      tabItem(tabName="About",
              box(
                solidHeader=TRUE,
                color='blue',
                width=20,
                title="Spatial Patterns and Transcriptomics using Ligand Activity (SPaTuLA)",
                "For our Bioengineering Capstone project we have created a spatial RNA-seq simulator. This project was for Northeastern University's Bioengineering Capstone and was created by 
                Nishanth Chinnadurai, Joseph Kern, Lucia Loosbrock, Magdalene Namwira, and Kaitlyn Ramesh.",
                br(),
                br(),
                "Our algorithm (SPaTuLA) is a spatial transcriptomics simulator that
enables bioengineers to input a cell signaling network to generate
synthetic spatial gene expression data. It integrates a hybrid AB-RD model
with the gene expression simulator, SymSim. We tested SPaTuLa on
theoretical and biological cell circuits to demonstrate its potential in
generating accurate tissue patterns.",
                br(),
                br(),
                "To take a look at our model code, head to our github page",
                br(),
                actionButton("github", "Github", onclick ="window.open('https://github.com/jtkern/ST-Capstone?tab=readme-ov-file', '_blank')"),
              )),
      
      #tutorial page
      tabItem(tabName='Tutorial',
              box(
                solidHeader=TRUE,
                status='primary',  
                width=6,
                title="Get started with SPaTuLA tutorial video",
                HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/DY7EMfWPvHo?si=3BsDHk5BigbzIx6e" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>'),
              ),
              box(
                solidHeader=TRUE,
                status='primary',
                width=6,
                title="Directions",
                "1. Navigate to SPaTuLA Simulator Input page",
                br(),
                "2. Enter desired inputs for SPaTuLA model (click 'Input Guide' for more details about each input)",
                br(),
                "3. Create a cell network",
                br(),
                "4. Input desired saving location and file names for gene expression matrix and animations",
                br(),
                "5. Enter any desired 'advanced' inputs (these inputs do not need to be changed but can be used to fine tune your network/simulator output)",
                br(),
                "6. Click 'Run' button to run the simulator",
                br(),
                "7. Navigate to SPaTuLA Simulator Output page to see graphical outputs of your simulation",
                br(),
                br(),
                actionButton("nav_inputs", "Input page")
              )),
      
      #page for model inputs
      tabItem(tabName="Inputs",
              
              fluidRow(
                box(
                  solidHeader=TRUE,
                  status="success",
                  collapsible=TRUE,
                  collapsed=TRUE,
                  title='Instructions',
                  "Use our prefilled input parameters to model a toggle switch circuit",
                  br(),
                  "Enter your own parameters to create your own simulation",
                  br(),
                  "Note that simulation outputs may take 2+ minutes to appear depending on complexity of network",
                  br(),
                  br(),
                  actionButton("switch_page",'Tutorial'),
                  actionButton("input_info", "Input Guide")
                ) 
              ),
              
              fluidRow(
                box(
                  solidHeader=TRUE,
                  status='primary', 
                  title="Input Parameters",
                  width=2,
                  numericInput("rho", "Cell density", "0.9"),
                  selectInput("cell_type_dist", "Choose cell-type distribution", choices = c("horizontal", "inwards radial","outwards radial", "multiradial", "vertical", "random", "uniform")),
                  selectInput("ligand_choice1", "Choose ligand 1 gradient", choices = c("inwards","outwards", "point", "uniform", "random", "angle", "multiradial")),
                  selectInput("ligand_choice2", "Choose ligand 2 gradient", choices = c("inwards","outwards", "point", "uniform", "random", "angle", "multiradial")),
                  numericInput("cell_dim", "Cell-cell matrix dimension", "30"),
                  numericInput("l1", "ligand threshold 1", "1"),
                  numericInput("l2", "ligand threshold 2", "20")
                ),
                
                box(
                  solidHeader=TRUE,
                  status='primary', 
                  title="Cell Network Input",
                  width=2,
                  #textInput("sender", "Enter sender cell", ""),
                  #textInput("receiver", "Enter receiver cell", ""),
                  numericInput("sender", "Enter sender cell", ""),
                  numericInput("receiver", "Enter receiver cell", ""),
                  selectInput("regulation2", "Choose regulation", choices = c(1,-1)),
                  numericInput("ligand", "Enter ligand", ""),
                  actionButton("add_data", "Add Data"),
                  actionButton("clear", "Clear Data"),
                ),
                
                box(
                  solidHeader=TRUE,
                  width=4,
                  title="Cell Interaction Input",
                  tableOutput("input_table")
                ),
                
                
                box(
                  solidHeader=TRUE,
                  status='primary', 
                  title="Save Output Files",
                  width=3,
                  "Mac users enter location as ~/.../.../",
                  br(),
                  "Windows users enter location as C:/.../.../",
                  br(),
                  textInput("file_location", "Enter desired file location", ""),
                  textInput("file_name1", "Enter gene expression matrix file name", "file1"),
                  textInput("file_name2", "Enter animation 1 name", "file2"),
                  textInput("file_name3", "Enter animation 2 name", "file3")
                ),
                
                box(  #run simulator
                  solidHeader=TRUE,
                  width=2,
                  status='success',
                  actionButton("run","Run Simulator")
                )
              ),
              
              fluidRow(
                box(
                  solidHeader=TRUE,
                  status='primary',  
                  width=6,
                  collapsible = TRUE,
                  collapsed = TRUE,
                  title="Advanced Inputs",
                  numericInput("g0", "Transcription factor basal expression rate","0.1"),
                  numericInput("g1", "Scaling factor for strength of regulation","5"),
                  numericInput("k", "Degradation rate of transcription factor", "0.1"),
                ), 
              )
      ),
      
      #page for model outputs
      tabItem(tabName="Outputs",
              
              box(
                solidHeader=TRUE,
                status='primary', 
                width=4,
                title='Gene Expression Heatmap',
                plotOutput("gene_exp_mat")
              ),
              
              box(
                solidHeader=TRUE,
                status='primary',
                width=4,
                title='Cell State Distribution',
                plotOutput("position_mat")
              ),
              
              box(
                solidHeader=TRUE,
                status='primary',
                width=4,
                title='PCA',
                plotOutput("PCA")
              )
              
      )
    )
  )
)

###################################

# Define server logic
server <- function(input, output, session) {
  
  
  observeEvent(input$switch_page, {    #switch to tutorial page on input page
    updateTabItems(session, "tabs", selected = "Tutorial")
  })
  
  observeEvent(input$nav_inputs, {    #switch to inputs page from tutorial page
    updateTabItems(session, "tabs", selected = "Inputs")
  })
  
  observeEvent(input$input_info, {     #input guide popup
    showModal(modalDialog(
      title = "Input Guide",
      "Sender cell: cell that is sending out the ligand (enter numeric)",
      br(),
      "Reciever cell: cell receiving input ligand (enter numeric)",
      br(),
      "Regulation: activation via ligand (1) or inhibition via ligand interaction (-1)",
      br(),
      "Ligand: ligand in reference to above inputs",
      br(),
      "Cell density: density of the cells in the system (matrix)",
      br(),
      "Cell-type distribution: how cells will be physically distributed across the matrix",
      br(),
      "Ligand gradient 1: direction of the gradient for ligand 1",
      br(),
      "Ligand gradient 2: direction of the gradient for ligand 2",
      br(),
      "Cell-cell matrix dimension: dimension of your cell matrix (i.e. enter 30 and will simulate a 30x30 matrix)",
      br(),
      "Ligand threshold 1: threshold for state change for ligand 1",
      br(),
      "Ligand threshold 2: threshold for state change for ligand 2",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  ### input matrix
  #input matrix that creates table using tibble function
  #change sender and reciever to character is using characters/text input
  input_table <- reactiveValues(
    table1 = tibble(sender = numeric(), receiver = numeric(), regulation = character(), ligand = numeric())
  )
  
  #add row to tibble
  observeEvent(input$add_data, {
    input_table$table1 <- input_table$table1 |> 
      add_row(sender = input$sender, receiver = input$receiver, regulation = input$regulation2, ligand = input$ligand)
  })
  
  #clear tibble
  observeEvent(input$clear, {
    input_table$table1 = tibble(sender = numeric(), receiver = numeric(), regulation = character(), ligand = numeric())
  })
  
  df_inputs <- as.data.frame(table1)  #save user input tibble as a data frame
  
  #outputs tibble to app
  output$input_table <- renderTable({
    input_table$table1
  })
  
  ##run SPaTuLA model
  observeEvent(input$run, {    #Run Simulator button
    
    updateTabItems(session, "tabs", selected = "Outputs")
    
    #define all input variables
    n=as.numeric(input$cell_dim)    # number grid points in each dimension, user defined
    rho = input$rho     #cell density, user defined
    mean_threshold = c(0.1, 1)
    dc = c(100, 5)
    k = c(0.1, 0.1)
    D = c(5,1)
    num_ligands = 2
    t_w = 1
    p = 1
    t_total = 150
    nframe = 151
    dt = 0.01
    param = c(g0=input$g0, g1=input$g1, k=input$k, n=5, l1=input$l1, l2=input$l2) #hill function, user input
    num_type = 2
    cell_type_choice=input$cell_type_dist
    centers = NULL
    props = c(0.5, 0.5)
    ligand_choice=c(input$ligand_choice1,input$ligand_choice2)
    center = c(n/2,n/2)
    conc_max = rep(10, 10)
    angle = NULL
    ngenes = 1000
    nevf = 100
    sim_gene_expr = FALSE
    saving_location= input$file_location
    file_name1=input$file_name1
    file_name2=input$file_name2
    file_name3=input$file_name3
    
    #use for now
    cell_net = data.frame(sender = c(1,2,1,2),
                          receiver=c(2,1,1,2),
                          regulation=c(-1,-1,1,1),
                          ligand=c(1,2,1,2))
    
    # cell_net = data.frame(sender = df_inputs[,1],           #This is gonna need troubleshooting
    #                       receiver = df_inputs[,2],
    #                       regulation = df_inputs[,3],
    #                       ligand = df_inputs[,4])
    
    ##################### SPaTuLA functions
    spatula_res = run_SPaTuLA(n=n, rho=rho, mean_threshold = mean_threshold, 
                              dc=dc, k=k, D=D,
                              num_ligands=num_ligands, t_w=t_w, p=p, 
                              t_total = t_total, dt=dt,
                              nframe=nframe, param=param, cell_net=cell_net, seed=0, 
                              num_type=num_type,
                              cell_type_choice=cell_type_choice, center=center, centers=centers, 
                              props=props, ligand_choice=ligand_choice, 
                              conc_max=conc_max, angle=angle, ngenes = ngenes, nevf = nevf, 
                              sim_gene_expr = sim_gene_expr)
    
    results = spatula_res[[1]]
    
    sim_grid = expand.grid(X=seq(n), Y=seq(n))
    frames = seq(nframe)
    l1_con_frames = cbind(rep(as.numeric(sim_grid$X), nframe),
                          rep(as.numeric(sim_grid$Y), nframe),
                          c(t(results$con_frames[[1]])), # need to transpose this matrix first!!
                          sort(rep(frames, n**2)))
    colnames(l1_con_frames) = c("X", "Y", "ligand1", "frame")
    
    l2_con_frames = cbind(rep(as.numeric(sim_grid$X), nframe),
                          rep(as.numeric(sim_grid$Y), nframe),
                          (c(t(results$con_frames[[2]]))), # need to transpose this matrix first!!
                          sort(rep(frames, n**2)))
    colnames(l2_con_frames) = c("X", "Y", "ligand2", "frame")
    
    tf_frames = cbind(rep(results$ind[,1], nframe),
                      rep(results$ind[,2], nframe),
                      c(t(results$tf_frames)),
                      sort(rep(frames, as.integer(n*n*rho))))
    colnames(tf_frames) = c("X", "Y", "tf", "frame")
    
    state_frames = cbind(rep(results$ind[,1], nframe), 
                         rep(results$ind[,2], nframe),
                         c(t(results$states)),
                         sort(rep(frames, as.integer(n*n*rho))))
    colnames(state_frames) = c("Xs", "Ys", "state", "frame")
    
    custom_colors <- c(  #could make user input 
      "1" = "#e0f2f1", # Light teal for Cell 1, Inactive
      "2" = "#80cbc4", # Medium teal for Cell 1, Active, Not Producing Ligand
      "3" = "#00695c", # Dark teal for Cell 1, Active, Producing Ligand
      "4" = "#ffe0b2", # Light peach for Cell 2, Inactive
      "5" = "#ffb74d", # Medium peach for Cell 2, Active, Not Producing Ligand
      "6" = "#e65100"  # Dark peach for Cell 2, Active, Producing Ligand
      # "7" = "#d9f0d3", # Light green for Cell 3, Inactive
      # "8" = "#a6d96a", # Medium green for Cell 3, Active, Not Producing Ligand
      # "9" = "#238443"  # Dark green for Cell 3, Active, Producing Ligand
    )
    
    g1 <- ggplot() + 
      # (1) Ligand concentration gradient
      geom_tile(data=as.data.frame(l1_con_frames), aes(X, Y, fill=ligand1)) + 
      scale_fill_viridis(option="B", name="[Ligand 1]") + # color scale for `con` variable
      
      # (2) Cell state distribution
      geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=as.factor(state)), size=2) + 
      scale_colour_manual(values = custom_colors, name="Cell State", limits = as.character(seq(num_type*3))) + 
      theme_bw() +
      theme(legend.position = "right") + # legend position
      transition_states(frame)
    # animate(g1, nframes=nframe)
    # animation <- animate(g1, nframes = nframe, renderer = file_renderer(saving_location, file_name1, overwrite = TRUE))
    # anim_save(paste0(saving_location, file_name1))
    
    g2 <- ggplot() + 
      # (1) Ligand concentration gradient
      geom_tile(data=as.data.frame(tf_frames), aes(X, Y, fill=tf)) + 
      scale_fill_viridis(option="A", name="[TF]") + 
      
      # (2) Cell state distribution
      geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=as.factor(state)), size=2) + 
      scale_colour_manual(values = custom_colors, name="Cell State", limits = as.character(seq(num_type*3))) + 
      theme_bw() +
      theme(legend.position = "right") + # legend position
      transition_states(frame)
    # animate(g2, nframes=nframe)
    # animation <- animate(g2, nframes = nframe, renderer = file_renderer(saving_location, file_name2, overwrite = TRUE))
    # anim_save(paste0(saving_location, file_name2))
    
    save(g1, file=paste0(saving_location, file_name1, ".rda"))
    save(g2, file=paste0(saving_location, file_name2, ".rda"))
    
    #PCA
    symsim_results = spatula_res[[2]]
    
    # Cell state distribution (final frame)
    state_frame = cbind(results$ind, results$states[nframe,])
    colnames(state_frame) = c("Xs", "Ys", "state")
    
    final_frame <- ggplot() + 
      # cell state distribution
      geom_point(data=as.data.frame(state_frame), aes(Xs, Ys, colour=as.factor(state)), size=2) + 
      theme(legend.position = "right", ) + 
      guides(colour=guide_legend(title="Cell State"))
    
    # Gene expression heatmaps
    gene_exprs_heatmap = DEG_analysis(results, symsim_results, nframe)
    
    ####################################################
    
    output$gene_exp_mat = renderPlot({
      gene_exprs_heatmap
    })
    
    output$position_mat = renderPlot({
      final_frame
    })
    
    output$PCA = renderPlot({
      spatula_res[[3]]
    })
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
