# RShiny app for UI and server; needs dynamic_helper.R function
# Changes to be made to adapt in this script have been commented with 'Provide'

# # Install + load libraries as required ----

if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
library(shiny)

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

if (!requireNamespace("bslib", quietly = TRUE)) {
  install.packages("bslib")
}
library(bslib)

if (!requireNamespace("shinythemes", quietly = TRUE)) {
  install.packages("shinythemes")
}
library(shinythemes)

if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  install.packages("shinyWidgets")
}
library(shinyWidgets)

if (!requireNamespace("decisionSupport", quietly = TRUE)) {
  install.packages("decisionSupport")
}
library(decisionSupport)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)  # For reading and writing CSV files

if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}
library(ggridges)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)
if (!requireNamespace("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}
library(ggtext)
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  install.packages("ggh4x")
}
library(ggh4x)
# Provide specific packages you have used other than the ones mentioned above
if (!requireNamespace("ggtext", quietly = TRUE)) {
  install.packages("ggtext")
}
library(ggtext)

if (!requireNamespace("png", quietly = TRUE)) {
  install.packages("png")
}
library(png)

if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}
library(grid)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

if (!requireNamespace("lubridate", quietly = TRUE)) {
  install.packages("lubridate")
}
library(lubridate)

# Load functins and inputs ----
# Provide Location of DA model script, dynamic-helper and funding-server scripts
# source("functions/saveLoad-module.R")
#source("functions/DA_for_exploring_funding_effects_data_visualisation.R")
source("functions/onion_model.R")
source("functions/dynamic-helper.R")
Source("functions/funding_server.R") # Provide this function when you have funding data entered in the excel sheet in data folder

# Provide Location of excel workbook containing the input parameters (prepared for the dynamic-helper)
file_path_vars <- "data/input_parameters.xlsx"
sheet_meta <- readxl::read_excel(file_path_vars, sheet = "sheet_names",
                                 col_types = c("text", "text"))
sheet_names <- sheet_meta$sheet_names
sheet_icons <- setNames(sheet_meta$icon, sheet_meta$sheet_names)


# Interface ----
ui <- fluidPage(
  
  theme = bs_theme(version = 5,
                   bootswatch = 'flatly',
                   base_font = font_google("Roboto")), 
  # Set actual browser tab title and favicon
  tags$head(
    tags$title("Agriculture Decision Support Tool"),
    tags$link(rel = "shortcut icon", href = "INRES.png"),
    
    tags$style(HTML("
    /* Scroll wrapper: scrolls horizontally *and* vertically only when needed */
    .scroll-xy {
      overflow-x: auto;                 /* left–right scroll  */
      overflow-y: auto;                 /* top–bottom scroll  */
      -webkit-overflow-scrolling: touch;/* smooth on iOS      */
      max-height: 80vh;                 /* optional: stop it taking more than
                                         80 % of the viewport height       */
  }
  
  /* Keep any Shiny plot inside that wrapper from shrinking */
  .scroll-xy .shiny-plot-output {
    min-width: 900px;                 /* choose your desktop width */
  }
                    ")
    )
  ),
  
  tags$div(
    style = "display:flex; align-items:center;justify-content:space-between;
      width: 100% !important; margin: 20px; padding: 0 15px;
      box-sizing: border-box; background-color: #f2f2f2;",
    
    # tags$a(href = "https://www.uni-bonn.de", target = "_blank",
    tags$img(src = "UniBonnHortiBonn_logo_transparent.png", height = "100px",
             style = "margin-left: auto; max-width: 20%; height: auto; cursor: pointer;"),
    # ),
    # Provide Title of the DA model
    tags$h2(tags$div("Climate Impact Forcasting:"),
            tags$div("conventional onion growing in 2075"),
            style = "text-align: center; flex-grow: 1;"),
    # Provide Project Logo
    # tags$a(href = "https://www.uni-bonn.de", target = "_blank",
    tags$img(src = "mlv-logo.png", height = "100px",
             style = "margin-right: auto; max-width: 30%; height: auto; cursor: pointer;")
    # ),
  ),
  
  
  ## Sidebar ----
  sidebarLayout(
    sidebarPanel(width = 4,
                 style = "height: 100%; overflow-y: auto",
                 
                 accordion(
                   id = "collapseSidebar",
                   open = FALSE,
                   
                   div(
                     class = "text-center",
                     actionButton("run_simulation", "Run Model",
                                  icon = icon("play"), class = "btn-primary")
                   ),
                   br(),
                   
                   ### Save/Load functionality ----
                   # saveLoadUI("savemod"),
                   accordion_panel(
                     title = "Save / Load project", icon = icon("folder-open"),
                     tagList(
                       textInput("state_name", "Project name"),
                       actionButton("save_btn",  label = tagList(icon("floppy-disk"),  "Save"  ), class = "btn btn-dark"),
                       
                       br(), br(),
                       selectInput("state_picker", "Saved versions", choices = NULL),
                       
                       fluidRow(
                         column(6, actionButton("load_btn",   tagList(icon("rotate"),  "Load"  ), class = "btn btn-secondary")),
                         column(6, actionButton("delete_btn", tagList(icon("trash"),   "Delete"), class = "btn btn-secondary"))
                       ),
                       hr(),
                       downloadButton("download_csv", label = tagList(icon("download"), "Download current inputs (.csv)"))
                     )
                   ),
                  
                   ### Expertise filter ----
                   accordion_panel(
                     title = "Expertise categories",
                     icon = icon("clipboard-question"),
                     tagList(
                       tags$h5(
                         "Expertise categories",
                         tags$span(
                           icon("circle-question"),
                           title = "Select your main expertise to view and edit only relevant variables.\nNot selecting any box shows all variables.\nDefaults apply to unselected categories in simulations.",
                           style = "cursor: help; margin-left: 8px;"
                         )
                       ),
                       uiOutput("category_filter_ui")
                     )
                   ),
                   
                   ### Funding scheme ----
                   accordion_panel(
                     title = "Funding schemes", icon = icon("euro-sign"),
                     create_funding_ui("funding")
                   ),
                   br(),
                   ### Dynamic elements ----
                   uiOutput("dynamic_element_ui")
                   
                 )
                 
    ),
    
## Main Panel ----
    mainPanel(width = 8,
              ###Tool description ----
              # Provide brief explanation of the DA model
              tags$h6(
                "This app simulates the present value of converting a treeless arable field into an alley cropping system with apple trees, based on a real farm in Germany.",
                tags$br(),
                tags$br(),
                "Use the tabs on the left to adjust variable ranges based on your local conditions or design goals.",
                tags$br(),
                tags$br(),
                "Click ‘Run model’ to perform a Monte Carlo simulation using random combinations from your defined ranges.You can save/load inputs, and once the model runs, results will appear below and you can save these figures.",
                tags$br(),
                tags$br(),
                "In the ‘Funding schemes’ tab, select any relevant funding options for your region.",
                tags$br(),
                #"DeFAF-suggested funding for German agroforestry: Annual support of 600 € per ha of wooded area and investment costs are to be funded at 100 % for first 10 ha of wooded area, 80 % for the next 10 ha, 50 % for additional area.",
                #tags$br(),
                "We welcome your feedback. For details on the onion model please feel free to contact",
                tags$a(href = "mailto:pauldaiber1996@gmail.com", "Paul Daiber"), ". For questions regarding the tool, you can reach out to", tags$a(href = "mailto:pkasargo@uni-bonn.de", "Prajna Kasargodu Anebagilu"), "or", tags$a(href = "mailto:afuelle1@uni-bonn.de", "Adrain Fuelle."),
              ),
              # tags$a(
              #   "Click here for latest info on Sustainable Farming Incentive",
              #   href = "https://www.gov.uk/government/publications/sustainable-farming-incentive-scheme-expanded-offer-for-2024",
              #   target="_blank",
              #   class = "my-btn"
              # ),
              br(), br(),
              
              #Debug:
              #titlePanel("Selected Financial Supports"),
              #create_funding_ui("funding"),        # UI part from the module
              #uiOutput("financial_support_links"),
              
              
              tags$h4("Selected Financial Supports"),
              tableOutput("summary"),
              
              #Debug:
              #verbatimTextOutput("summary"), 
              # uiOutput("funding-financial-support"),
              
              ### Plots ----
              br(), br(),
              div(class = "scroll-xy",
                  plotOutput("plot1_ui", height = "550px"),
              ),
              br(),
              uiOutput("plot1_dl_ui"),
              br(), br(),br(), br(),
              
              div(class = "scroll-xy",
                  plotOutput("plot2_ui", height = "550px"),
              ),
              br(),
              uiOutput("plot2_dl_ui"),
              br(), br(),br(), br(),
              
              div(class = "scroll-xy",
                  plotOutput("plot3_ui", height = "700px"),
              ),
              br(),
              uiOutput("plot3_dl_ui"),
              br(), br(),br(), br(),
              
              div(class = "scroll-xy",
                  plotOutput("plot4_ui", height = "550px"),
              ),
              br(),
              uiOutput("plot4_dl_ui"),
              br(), br(),br(), br(),
              
              div(class = "scroll-xy",
                  plotOutput("plot5_ui", height = "550px"),
              ),
              br(),
              uiOutput("plot5_dl_ui"),
              br(), br(),br(), br(),
              
              # div(class = "scroll-xy",
              # plotOutput("plot6_ui", height = "550px"),
              # ),
              # br(),
              # uiOutput("plot6_dl_ui"),
              # br(), br(),br(), br(),
              # 
              # div(class = "scroll-xy", plotOutput("plot7_ui", height = "550px"),),
              # br(),
              # uiOutput("plot7_dl_ui"),
              # br(), br(),br(), br(),
              
              div(class = "scroll-xy", 
                  plotOutput("plot8_ui", height = "550px"),
              ),
              br(),
              uiOutput("plot8_dl_ui"),
              
              ### Requisites ----
              # Provide Funding declaration
              tags$img(src = "Funding_declaration.png", height = "100px",
                       style = "margin-right: auto; max-width: 100%; height: auto; cursor: pointer;"),
              #Provide tool usage disclaimer
              tags$p(
                  tags$a("Disclaimer", href = "https://agroreforest.eu/reforest-tools-disclaimer/",
                target = "_blank"),
                     " | ",  
              #Provide the correct link once the app and codes are hosted in HortiBonn repo
              tags$a("View Source", href = "https://github.com/PrajnaKARai/Onion_DA", # temp link
                target = "_blank")
                  ),
              br(), br(),br(), br(),
              
    )
  )
  
)


# Server ----
server <- function(input, output, session) {
  
  ## Dynamic funding module ----
  funding <- funding_server("funding")   # returns a list of reactives
  
  output$`funding-financial-support` <- renderUI({
    funding$financial_support_links
  })
  # output$summary <- renderPrint({
  #   result$category_totals()          
  output$summary <- renderTable({
    
    # Get the full funding totals (gov + private) as a named list
    total_funding <- funding$total_funding_with_private()
    
    # Debug: data frame for table output ### remove for the final or can be displayed in the mainPanel too - upto @Adrain
    data.frame(
      `Funding Category` = str_to_title(str_replace_all(str_remove(names(total_funding), "_c$"), "_", " ")),
      `Total Financial Support` = round(unname(total_funding), 2),
      check.names = FALSE,
      row.names = NULL
    )
  })
  ## Helper for safe extraction from named vector
  safe_get <- function(vec, name) {
    if (is.null(vec) || length(vec) == 0 || is.na(vec[name])) return(0)
    if (! name %in% names(vec)) return(0)
    as.numeric(vec[name])
  }
  
  ## Dynamic expertise-filter module ----
  # helper that sanitises category names into safe IDs
  sanitize <- function(x) gsub("[^A-Za-z0-9]", "_", x)
  
  # all categories across every sheet
  categories <- reactive({
    cats <- unique(unlist(lapply(excelData(), function(df) df$Expertise)))
    cats <- cats[!is.na(cats) & cats != ""]
    trimws(unique(unlist(strsplit(cats, ";"))))
  })
  
  # filter bar UI
  output$category_filter_ui <- renderUI({
    if (length(categories()) == 0) return(NULL)
    tagList(
      lapply(categories(), function(cat){
        checkboxInput(
          paste0("cat_", sanitize_id(cat)), cat, value = FALSE)
      })
    )
  })
  
  ## Dynamic UI inputs ----
  
  # read in input xlsx file
  excelData <- reactive({
    sheet_number <- seq_along(sheet_names)+1
    all_sheets <- lapply(sheet_number, function(sht) {
      readxl::read_excel(file_path_vars, sheet = sht,
                         col_types = c("text", "numeric", "numeric", "text", "text", "text", "text", "guess", "guess", "text", "text")
      )
    })
    names(all_sheets) <- sheet_names
    all_sheets
  })
  
  
  # util: turns a category vector into a JS condition 
  ### render but hide unchecked expertise categories - default show-all ----
  panel_condition <- function(cat_vec) {
    cat_vec <- trimws(cat_vec)
    cat_vec <- cat_vec[cat_vec != "" & !is.na(cat_vec)]
    if (length(cat_vec) == 0) return("true")
    
    cat_ids <- sprintf("input['cat_%s']", sanitize_id(cat_vec))
    
    cat_show_all <- paste0(
      "Object.keys(input).filter(k => k.startsWith('cat_')).",
      "every(k => input[k] === false)"
    )
    
    sprintf("(%s) || (%s)",            # show when *no* cat box ticked
            cat_show_all,              #…or any matching cat ticked
            paste(cat_ids, collapse = ' || '))
  }
  
  
  output$dynamic_element_ui <- renderUI({
    
    data_list   <- excelData()
    sheet_names <- names(data_list)
    
    # build one accordion panel per sheet
    # the elements are generated via the external function create_ui_element()
    panels <- lapply(seq_along(data_list), function(j) {
      
      sheet <- data_list[[j]]
      
      cats  <- unique(trimws(unlist(strsplit(sheet$Expertise %||% "", ";|,"))))
      cats  <- cats[cats != ""]
      
      ui_elems <- lapply(seq_len(nrow(sheet)), function(i) {
        create_ui_element(sheet[i, ])
      })
      
      conditionalPanel(
        condition = panel_condition(cats),   # hide panel is empty
        accordion_panel(
          title = sheet_names[j],
          icon  = icon(sheet_icons[[ sheet_names[j] ]] %||% "circle-dot"),
          tagList(ui_elems)
        )
      )
    })
    
    tagList(panels)   # render the list
  })
  
  
  ## Save, Load and Delete module
  all_inputs <- reactive({
    names(input)[grepl("(_c$|_p$|_t$|_n$|_cond$)", names(input))]
  })
  
  current_input_table <- reactive({
    variables <- all_inputs()
    # # source("functions/Walnut_grain_veg_tub_mcsim-only.R", local = T)
    # 
    # message("Accessing user input estimates from the interface...")
    # 
    # # 1. Gather current widget values
    # exclude_inputs <- c("collapseSidebar", "save", "load", "delete",
    #                     "confirm_delete", "admin_selected_user",
    #                     "project_name", "version_select", "delete_version_select")
    # 
    # # variables <- setdiff(
    #   names(input)[grepl("(_c$|_p$|_t$|_n$|_cond$)", names(input))],
    #   exclude_inputs
    # )
    
    lower_values <- sapply(variables, function(v) {
      val <- input[[v]]
      if (length(val) == 1) as.numeric(val) else as.numeric(val[1])
    })
    upper_values <- sapply(variables, function(v) {
      val <- input[[v]]
      if (length(val) == 1) as.numeric(val) else as.numeric(val[2])
    })
    
    # 2. Re-read Excel (keeps original bounds & distributions)
    all_sheets <- excelData()            # list of data-frames
    input_file <- bind_rows(all_sheets)  # one big table
    
    # Overwrite lower/upper with current UI inputs
    input_file <- input_file %>%
      left_join(
        tibble(variable = variables,
               lower    = lower_values,
               upper    = upper_values),
        by = "variable",
        suffix = c("", ".new")
      ) %>%
      mutate(
        lower = coalesce(lower.new, lower),
        upper = coalesce(upper.new, upper)
      ) %>%
      select(-ends_with(".new"))
    
    # 3. Append funding scalars
    # View(input_file)
    #print(1)
    
    funding_names <- 
      c("funding_onetime_percentage_initial_cost_schemes_c", "annual_funding_schemes_c",
        "funding_onetime_percentage_consult_schemes_c","funding_onetime_per_tree_schemes_c",
        "funding_onetime_per_m_treerow_schemes_c", "funding_onetime_per_m_hedgerow_schemes_c","annual_funding_per_m_schemes_c",
        "annual_funding_per_tree_schemes_c", "funding_onetime_schemes_c",
        "onetime_external_percentage_incost_schemes_c","onetime_external_percentage_consult_schemes_c",
        "funding_onetime_per_ha_schemes_c", "onetime_external_support_c", "annual_external_support_c")
    funding_df <- data.frame(variable = funding_names,
                             lower = 0,
                             upper = 0,
                             distribution = "const")
    
    try(total_funding <- funding$total_funding_with_private())
    
    if ("total_funding" %in% ls()) {
      #print(2)
      
      input_file <- 
        data.frame(variable = names(total_funding),
                   lower = unname(total_funding),
                   upper = unname(total_funding),
                   distribution = "const") %>% 
        bind_rows(input_file, .)
      
      remain <- funding_names[!(funding_names %in% input_file$variable)]
      input_file <- funding_df %>% 
        filter(variable %in% remain) %>% 
        bind_rows(input_file, .)
      
      # View(input_file)
    }else {
      input_file <- bind_rows(input_file, funding_df)
    }
    
    #View(input_file)
    
    # # 4. Save UI snapshot (optional)
    # saveRDS(list(sheet_names, input_file), "data/Walnut_grain_veg_tub_ui_updated.RDS")
    
    # 5. clean-up: keep only numeric rows
    input_file <- input_file %>%
      filter(
        !is.na(lower), !is.na(upper),
        is.finite(lower), is.finite(upper)
      )
    
    # write.csv(input_file,"data/input_table.csv",row.names = F)
    
    input_file
  })
  
  ## Save/Load functionality ----
  # saveLoadServer("savemod", current_input_table)
  # Provide Folder name instead of the current 'Germany' to store user saves
  get_base_dir <- function() {
    if (Sys.info()[["sysname"]] == "Windows")
      "user-states/Germany_onion"
    else
      "/srv/shiny-app-data/user-states/Germany_onion"
  }
  
  get_user_dir <- function() {
    uid <- session$user
    safe_uid <- if (is.null(uid) || uid == "") "anon"
    else gsub("[^A-Za-z0-9_.-]", "_", uid)
    dir <- file.path(get_base_dir(), safe_uid)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir
  }
  
  timestamp_name <- function(raw) {
    paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),"_",
           gsub("[^A-Za-z0-9_.-]", "_", raw), ".rds")
  }
  
  observeEvent(input$save_btn, {
    dir  <- get_user_dir()
    files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
    
    if (length(files) >= 5) {
      showModal(modalDialog("You already have five versions. Delete one first.",
                            easyClose = TRUE))
      return()
    }
    
    req(nzchar(input$state_name))
    saveRDS(
      list(input_table = current_input_table(),
           raw_inputs  = reactiveValuesToList(input)),
      file.path(dir, timestamp_name(input$state_name))
    )
  })
  
  saved_files <- reactiveFileReader(
    2000, session, get_user_dir(),
    function(dir) sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE),decreasing = T)
  )
  
  observe({
    updateSelectInput(session, "state_picker",
                      choices = basename(saved_files()))
  })
  
  observeEvent(input$load_btn, {
    req(input$state_picker)
    obj <- readRDS(file.path(get_user_dir(), input$state_picker))
    bslib::accordion_panel_open("collapseSidebar",TRUE,session)
    vals <- obj$raw_inputs
    
    restore_one <- function(id, val) {
      if (is.null(val)) return()
      switch(class(val)[1],
             numeric   = updateNumericInput(session, id, value = val),
             integer   = updateNumericInput(session, id, value = val),
             character = updateTextInput   (session, id, value = val),
             logical   = updateCheckboxInput(session, id, value = val),
             factor    = updateSelectInput (session, id, selected = as.character(val)),
             # length-2 numeric == slider
             { if (is.numeric(val) && length(val) == 2)
               updateSliderInput(session, id, value = val) }
      )
    }
    
    # ordinary widgets
    lapply(names(vals), \(id) try(restore_one(id, vals[[id]]), silent = TRUE))
    
    # funding module widgets  (country + state first, the rest after rebuild)
    ## doesn't work cleanly yet - load button needs to be pressed twice
    ns <- NS("funding")   # helper to prepend "funding-"
    
    # (a) push country and state immediately 
    try(updateSelectInput(session, ns("country"),
                          selected = vals[[ns("country")]]), silent = TRUE)
    try(updateSelectInput(session, ns("state"),
                          selected = vals[[ns("state")]]),   silent = TRUE)
    
    # (b) *once* the state really is set, restore the rest
    observeEvent(input[[ns("state")]], {
      if (!identical(input[[ns("state")]], vals[[ns("state")]])) return()
      
      try(updateSelectInput(session, ns("one_schemes"),
                            selected = vals[[ns("one_schemes")]]), silent = TRUE)
      try(updateSelectInput(session, ns("annual_schemes"),
                            selected = vals[[ns("annual_schemes")]]), silent = TRUE)
      try(updateNumericInput(session, ns("onetime_private"),
                             value = vals[[ns("onetime_private")]]),  silent = TRUE)
      try(updateNumericInput(session, ns("annual_private"),
                             value = vals[[ns("annual_private")]]),   silent = TRUE)
    }, once = TRUE, ignoreInit = FALSE)
  })
  
  observeEvent(input$delete_btn, {
    req(input$state_picker)
    unlink(file.path(get_user_dir(), input$state_picker))
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("current_input_", Sys.Date(), ".csv"),
    content  = function(file) write_csv(current_input_table(), file)
  )
  
  
  
  ## Monte Carlo Simulation ----
  mcSimulation_results <- eventReactive(input$run_simulation, {
    input_file <- current_input_table()
    
    # 6. Run Monte-Carlo
    # Provide model_function
    decisionSupport::mcSimulation(
      estimate          = decisionSupport::as.estimate(input_file),
      model_function    = AF_benefit_with_Risks,
      numberOfModelRuns = input$num_simulations_c,
      functionSyntax    = "plainNames"
    )
    
  })
  
  ## Generating plots ----
  # helper to add title subtile caption etc
  add_meta <- function(p, title, subtitle = NULL, caption = NULL,
                       legend = "bottom") {
    
    p +
      labs(title = title, subtitle = subtitle, caption = caption) +
      theme(
        plot.title = element_textbox_simple(
          size   = 24,
          face   = "bold",
          width  = unit(1, "npc"),  # full plot width
          halign = 0.5,              # centered
          margin = margin(t = 6,b = 20)
        ),
        plot.subtitle = element_textbox_simple(
          size   = 18,
          width  = unit(1, "npc"),
          halign = 0.5,
          margin = margin(t = 6,b = 20)
      ),
    plot.caption  = element_textbox_simple(
      size   = 16,
      width  = unit(0.98, "npc"),
      halign = 0,              # left-aligned
      margin = margin(t = 6,b = 20),
      hjust = 0,
      vjust = 1
    ),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text     = element_text(size = 14, hjust = 0.5),
    legend.position = legend,
    plot.margin = margin(t = 50, r = 10, b = 50, l = 10, unit = "pt")
    
    )  }
  
  # download helper
  make_download <- function(id, plot_obj, filename, width = 13, height = 5, dpi = 300, scale = 2) {
    output[[id]] <- downloadHandler(
      filename = function() filename,
      content  = function(file) {
        # device is inferred from file extension; here "png"
        ggsave(file, plot_obj, width = width, height = height, dpi = dpi, scale = scale)
      }
    )
  }
  
  
  
  observeEvent(mcSimulation_results(), {
    mc_data <- mcSimulation_results()
    # Provide correct variables for plots
    plot1 <-
      decisionSupport::plot_distributions(mcSimulation_object = mc_data,
                                          vars = c("NPV_Treeless_System", "NPV_Agroforestry_System"),
                                          method = "boxplot",
                                          #method = "smooth_simple_overlay",
                                          old_names = c("NPV_Treeless_System", "NPV_Agroforestry_System"),
                                          new_names = c("Monoculture (baseline)", "Agroforestry with current funding"),
                                          x_axis_name = "NPV (€)",
                                          y_axis_name = "Decision Options") |>
      add_meta(
        title    = "Figure 1. Probabilistic distributions of Net Present Value",
        subtitle = "Agroforestry intervention with current funding vs. conventional farming",
        caption  = "Figure 1 shows the comparison of Net Present Value (NPV) outcomes for agroforestry (Apple alley cropping) vs monoculture system (baseline). The x-axis displays NPV values (i.e., the sum of discounted annual cash flows). The higher and wider the box, the greater the potential return and variability in outcomes under that system."
      )
    
    plot2 <- decisionSupport::plot_distributions(
      mc_data, "NPV_decis_AF_ES3",
      method     = "smooth_simple_overlay",
      old_names  = "NPV_decis_AF_ES3",
      new_names  = "Agroforestry – Treeless",
      x_axis_name= "NPV (€)",
      y_axis_name= "Probability") |>
      add_meta(
        title    = "Figure 2. Distribution of the *incremental* NPV",
        subtitle = "Difference between agroforestry and treeless farming under identical conditions",
        caption  = "Figure 2 shows the NPV distributions of the decision to establish the apple alley cropping system
                as compared to the decision to continue with monoculture for the specified time (i.e., NPV agroforestry - NPV monoculture under identical conditions).
                The x-axis displays NPV values (i.e., the sum of discounted annual cash flows) and y-axis displays the probability of each NPV amount to occur (i.e., higer y-values indicate higher probability)"
        , legend = "none")
    
    plot3 <- decisionSupport::plot_distributions(
      mc_data,
      vars      = c("NPV_decis_no_fund", "NPV_decis_AF_ES3", "NPV_decis_DeFAF"),
      method    = "boxplot",
      old_names = c("NPV_decis_no_fund", "NPV_decis_AF_ES3", "NPV_decis_DeFAF"),
      new_names = c("Agroforestry without\nfunding - Treeless", "Agroforestry with\ncurrent funding - Treeless", "Agroforestry with\nDeFAF-suggested funding - Treeless"),
      x_axis_name = "NPV (€)",
      y_axis_name = "Funding Options") |>
      add_meta(
        title    = "Figure 3. Net Present Value (NPV) Outcomes Across Funding Schemes for Apple Alley Cropping",
        subtitle = "Agroforestry intervention with, without and DeFAF-suggested funding",
        caption  = 'Figure 3 shows the comparison of net present value (NPV) outcomes for the decision of different agroforestry funding schemes. The x-axis displays NPV values (i.e., the sum of discounted annual cash flows); each colored boxplot represents a funding scheme, showing the range and distribution of simulation results from the probabilistic model.
        The higher and wider the box, the greater the potential return and variability in outcomes under that funding.
Scenarios involving funding (like DeFAF-suggested or EcoScheme3 and regional) generally show higher NPV ranges than the No funding, however it not necessarily better suggesting the current financial support is insufficient to sustain agroforestry.'
      )
    
    plot4 <- decisionSupport::plot_cashflow(
      mc_data, "AF_CF",
      x_axis_name = "",
      y_axis_name = "Annual cash-flow from Agroforestry (€)",
      color_25_75 = "navajowhite",
      color_5_95 = "green4",
      color_median = "darkblue",
      facet_labels = "") |>
      add_meta(
        title   = "Figure 4. Annual cash-flow of the agroforestry intervention", 
        subtitle = "Projected yearly cash-flow variability for an agroforestry system over time",
        caption = 'Figure 4 shows how annual cash-flow from an agroforestry intervention is expected to evolve, based on a probabilistic simulation. The shaded areas represent uncertainty ranges (from lower to upper quantiles), while the blue line shows the median outcome (expressed in €). While early years may involve negative cash flow, profitability tends to improve over time, with increasing stability. The graph highlights the long-term financial potential and risk spread of adopting agroforestry practices.'
      )
    
    plot5 <- decisionSupport::plot_cashflow(
      mc_data, "AF_CCF_ES3",
      x_axis_name = "",
      y_axis_name = "Cumulative cash-flow from Agroforestry (€)",
      color_25_75 = "navajowhite",
      color_5_95 = "green4",
      color_median = "darkblue",
      facet_labels = "") |>
      add_meta(
        title   = "Figure 5. Cumulative cash-flow of the agroforestry intervention", 
        subtitle = "Long-term cumulative cash-flow projection for an agroforestry system",
        caption = "Figure 5  illustrates how total cash-flow (expressed in €) accumulates over time from an agroforestry intervention, based on a range of simulated outcomes. The shaded areas represent uncertainty (spread of possible results), and the blue line indicates the median trajectory. Cumulative returns grow steadily over time, showing the long-term profitability potential of agroforestry. Despite initial variability, the system trends positively, reinforcing the case for agroforestry as a viable financial investment over the long run."
      )
    
    # plot6 <- decisionSupport::plot_cashflow(
    #   mc_data, "Cashflow_AF1_decision",
    #   x_axis_name = "",
    #   y_axis_name = "Annual cash-flow (€)",
    #   facet_labels = "") |>
    #   add_meta(
    #     title   = "Figure 6. Incremental annual cash-flow",
    #     subtitle= "Agroforestry minus baseline farming",
    #     caption = "Figure 6 shows the difference (expressed in €) between the annual balance of alley-cropping and continue farming without planting trees under identical real-world scenarios.")
    # 
    # plot7 <- decisionSupport::plot_cashflow(
    #   mc_data, "Cum_Cashflow_AF1_decision",
    #   x_axis_name = "",
    #   y_axis_name = "Cumulative cash-flow (€)",
    #   facet_labels = "") |>
    #   add_meta(
    #     title   = "Figure 7. Incremental cumulative cash-flow",
    #     subtitle= "Agroforestry minus baseline farming",
    #     caption = 'Figure 7 shows the cumulative difference (expressed in €) between the annual balance of alley-cropping and continue farming without planting trees under identical real-world scenarios.')
    # 
    
    # Send plots to UI
    output$plot1_ui <- renderPlot({ plot1 })
    make_download("download_plot1", plot1, "Figure1_NPV.png")
    output$plot1_dl_ui <- renderUI({
      downloadButton("download_plot1", "Download Figure 1")
    })
    
    output$plot2_ui <- renderPlot({ plot2 })
    make_download("download_plot2", plot2, "Figure2_Decision_NPV.png")
    output$plot2_dl_ui <- renderUI({
      downloadButton("download_plot2", "Download Figure 2")
    })
    
    output$plot3_ui <- renderPlot({ plot3 })
    make_download("download_plot3", plot3, "Figure3_Funding_NPVs.png")
    output$plot3_dl_ui <- renderUI({
      downloadButton("download_plot3", "Download Figure 3")
    })
    
    output$plot4_ui <- renderPlot({ plot4 })
    make_download("download_plot4", plot4, "Figure4_Annual_Cashflow.png")
    output$plot4_dl_ui <- renderUI({
      downloadButton("download_plot4", "Download Figure 4")
    })
    
    output$plot5_ui <- renderPlot({ plot5 })
    make_download("download_plot5", plot5, "Figure5_Cumulative_Cashflow.png")
    output$plot5_dl_ui <- renderUI({
      downloadButton("download_plot5", "Download Figure 5")
    })
    
    # output$plot6_ui <- renderPlot({ plot6 })
    # make_download("download_plot6", plot6, "Figure6_Incremental_Annual_CF.png")
    # output$plot6_dl_ui <- renderUI({
    #   downloadButton("download_plot6", "Download Figure 6")
    # })
    # 
    # output$plot7_ui <- renderPlot({ plot7 })
    # make_download("download_plot7", plot7, "Figure7_Incremental_Cumulative_CF.png")
    # output$plot7_dl_ui <- renderUI({
    #   downloadButton("download_plot7", "Download Figure 7")
    # })
    
    
    # Ask user whether to run EVPI (takes time!)
    showModal(modalDialog(
      title = "Run EVPI analysis?",
      "Do you want to assess the Expected Value of Perfect Information (EVPI)?
      This step may take a while, but you can explore the other graphs while the EVPI is processed.
      The EVPI graph will appear at the bottom of the page, below the last graph.",
      footer = tagList(
        modalButton("No"),
        actionButton("confirm_evpi", "Yes, run EVPI")
      )
    ))
    
    # Handle user confirmation to run EVPI
    observeEvent(input$confirm_evpi, {
      
      removeModal()  # remove popup
      
      # Try running EVPI only if it can return meaningful values
      tryCatch({
        evpi_input <- as.data.frame(cbind(
          mc_data$x,
          #Provide correct variable
          NPV_decision_AF1 = mc_data$y$NPV_decis_AF_ES3
        ))
        # Provide the NPV_decision variable to calculate EVPI
        evpi_result <- decisionSupport::multi_EVPI(evpi_input, "NPV_decis_AF_ES3")
        
        # saveRDS(evpi_input, "evpi_input_test.rds")
        # evpi_input <- readRDS("evpi_input_test.rds")
        
        # saveRDS(evpi_result, "evpi_result_test.rds")
        # evpi_result <- readRDS("evpi_result_test.rds")
        
        var_lookup <- bind_rows(excelData()) %>%
          filter(!is.na(variable), !is.na(name)) %>%
          distinct(variable, name) %>%
          deframe()
        #Provide correct variable
        plot8 <- plot_evpi(evpi_result, decision_vars = "NPV_decis_AF_ES3",
                           new_names = "") +
          scale_y_discrete(labels = var_lookup)
        
        plot8 <- plot8 |>
          add_meta(title = "Figure 8. EVPI for Each Variable",
                   subtitle = "Maximum amount worth paying for perfect information on each variable."
          )
        
        output$plot8_ui <- renderPlot({ plot8 })
        
        make_download("download_plot8", plot8, "Figure8_EVPI.png")
        
        output$plot8_dl_ui <- renderUI({
          downloadButton("download_plot8", "Download Figure 8")
        })
        
      }, error = function(e) {
        warning("EVPI plot skipped due to error: ", e$message)
        output$plot8_ui <- renderPlot({
          plot.new()
          text(0.5, 0.5, "There are no variables with a positive EVPI.\nGetting better information will \nnot reduce the level of uncertainty of the decision.", cex = 1.2)
        })
      })
    })
    
  })
  
}

shinyApp(ui = ui, server = server)
