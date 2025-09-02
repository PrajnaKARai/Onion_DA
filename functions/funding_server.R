# dynamic function to create funding UI based on the structure defined in funding_data.xlsx
# function to be called in app.R 
# Provide excel containing funding data as structured in funding_data1.xlsx
create_funding_ui <- function(id, data_path = "data/funding_data1.xlsx", sheet = 1) {
  ns <- NS(id) # use namespacing instead of concatenating values - ask @Adrain if this is okay...?####
  funding_data <- readxl::read_excel(data_path, sheet = sheet)
  
  tagList(
    selectInput(ns("country"), "Select Country", choices = unique(funding_data$Country)),
    uiOutput(ns("state_ui")),
    # Government aids/subsidies/ financial support
    
    uiOutput(ns("funding_one_ui")),      # One-time funding dropdown
    uiOutput(ns("funding_annual_ui")),   # Annual funding dropdown
    # private financial support from excel only
    uiOutput(ns("funding_private_onetime_ui")),  # One-time Private dropdown + numeric
    uiOutput(ns("funding_private_annual_ui")),    # Annual Private dropdown + numeric
    # private financial support diretc user entry
    uiOutput(ns("donation_private_onetime_ui")),  # One-time Private dropdown + numeric
    uiOutput(ns("donation_private_annual_ui")),    # Annual Private dropdown + numeric
    
    uiOutput(ns("financial_support_links")) # Link to the scheme
  )
}


funding_server <- function(id, data_path = "data/funding_data1.xlsx", sheet = 1) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    funding_data <- readxl::read_excel(data_path, sheet = sheet)
    
    # Normalize Funding Type column
    funding_data$Funding_Type_Clean <- tolower(trimws(funding_data$`Funding Type`))
    
    # Country-based state filtering
    filtered_states <- reactive({
      req(input$country)
      funding_data %>%
        filter(Country == input$country) %>%
        pull(State) %>%
        unique() %>%
        sort()
    })
    
    output$state_ui <- renderUI({
      selectInput(ns("state"), "Select State", choices = filtered_states())
    })
    
##############################################
## INSIDE moduleServer(...) ##################

# remember previously-selected items for each drop-down
prev <- reactiveValues(
  one_gov    = character(0),
  annual_gov = character(0),
  one_priv   = character(0),
  annual_priv= character(0)
)

# helper that pops up a modal if txt (vector) is non-empty
show_if_has_condition <- function(txt) {
  txt <- txt %>% na.omit() %>% trimws() %>% .[. != ""]
  if (length(txt) > 0) {
    showModal(modalDialog(
      title  = "Requirement(s) for this support",
      paste(txt, collapse = "\n\n"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  }
}

# ------- ONE-TIME GOV ------------------------------------
observeEvent(input$one_schemes, {
  new <- setdiff(input$one_schemes, prev$one_gov)   # what was just added?
  prev$one_gov <- input$one_schemes                 # update memory
  if (length(new)) {
    cond <- filtered_data() %>%
      filter(`UI display` %in% new) %>% pull(Condition)
    show_if_has_condition(cond)
  }
})

# ------- ANNUAL GOV --------------------------------------
observeEvent(input$annual_schemes, {
  new <- setdiff(input$annual_schemes, prev$annual_gov)
  prev$annual_gov <- input$annual_schemes
  if (length(new)) {
    cond <- filtered_data() %>%
      filter(`UI display` %in% new) %>% pull(Condition)
    show_if_has_condition(cond)
  }
})

# ------- ONE-TIME PRIVATE --------------------------------
observeEvent(input$one_private_schemes, {
  new <- setdiff(input$one_private_schemes, prev$one_priv)
  prev$one_priv <- input$one_private_schemes
  if (length(new)) {
    cond <- filtered_data() %>%
      filter(`UI display` %in% new) %>% pull(Condition)
    show_if_has_condition(cond)
  }
})

# ------- ANNUAL PRIVATE ----------------------------------
observeEvent(input$annual_private_schemes, {
  new <- setdiff(input$annual_private_schemes, prev$annual_priv)
  prev$annual_priv <- input$annual_private_schemes
  if (length(new)) {
    cond <- filtered_data() %>%
      filter(`UI display` %in% new) %>% pull(Condition)
    show_if_has_condition(cond)
  }
})

# -------- OPTIONAL: reset memory when country/state changes
observeEvent(c(input$country, input$state), {
  prev$one_gov    <- character(0)
  prev$annual_gov <- character(0)
  prev$one_priv   <- character(0)
  prev$annual_priv<- character(0)
})
##############################################

    
    
    
    # Data filtered by country and state (including schemes marked as 'Any')
    filtered_data <- reactive({
      req(input$country, input$state)
      funding_data %>%
        filter(Country == input$country) %>%
        filter(State == input$state | State == "Any")
    })
    
    # One-time GOV funding
    output$funding_one_ui <- renderUI({
      one_time_gov <- filtered_data() %>%
        filter(grepl("^one-time", Funding_Type_Clean) & !grepl("private", Funding_Type_Clean))
      
      selectInput(
        ns("one_schemes"),
        "Select Government One-time Funding Schemes",
        choices = one_time_gov$`UI display`,
        multiple = TRUE
      )
    })
    
    # Annual GOV funding
    output$funding_annual_ui <- renderUI({
      annual_gov <- filtered_data() %>%
        filter(grepl("^annual", Funding_Type_Clean) & !grepl("private", Funding_Type_Clean))
      
      selectInput(
        ns("annual_schemes"),
        "Select Government Annual Funding Schemes",
        choices = annual_gov$`UI display`,
        multiple = TRUE
      )
    })
    
    # One-time PRIVATE funding
    output$funding_private_onetime_ui <- renderUI({
      one_time_private <- filtered_data() %>%
        filter(grepl("^one-time", Funding_Type_Clean) & grepl("private", Funding_Type_Clean))
      # << hide entire widget when no private options are available >>
      if (nrow(one_time_private) == 0) return(NULL)
      
      tagList(
        selectInput(
          ns("one_private_schemes"),
          "Select One-time Private Support",
          choices = one_time_private$`UI display`,
          multiple = TRUE
        ),
        #numericInput(ns("onetime_private_custom"), "One-time funding from Crowdfunding, donations etc.", min = 0, max = 100000, value = 0)
      )
    })
    
    # Annual PRIVATE funding
    output$funding_private_annual_ui <- renderUI({
      annual_private <- filtered_data() %>%
        filter(grepl("^annual", Funding_Type_Clean) & grepl("private", Funding_Type_Clean))
    
    # << hide entire widget when no private options are available >>
    if (nrow(annual_private) == 0) return(NULL)
      
      tagList(
        selectInput(
          ns("annual_private_schemes"),
          "Select Annual Private Support",
          choices = annual_private$`UI display`,
          multiple = TRUE
        ),
        #numericInput(ns("annual_private_custom"), "Annual private donations received [money/ha]", min = 0, max = 5000, value = 0)
      )
    })
    
    
    # One-time PRIVATE funding
    output$donation_private_onetime_ui <- renderUI({
      # d_one_time_private <- filtered_data() %>%
      #   filter(grepl("^one-time", Funding_Type_Clean) & grepl("private", Funding_Type_Clean))

      tagList(
        numericInput(ns("onetime_private_custom"), "One-time funding from Crowdfunding, donations etc.", min = 0, max = 100000, value = 0)
      )
    })
    
    # Annual PRIVATE funding
    output$donation_private_annual_ui <- renderUI({
      # d_annual_private <- filtered_data() %>%
      #   filter(grepl("^annual", Funding_Type_Clean) & grepl("private", Funding_Type_Clean))
      
      tagList(
        numericInput(ns("annual_private_custom"), "Annual private donations received [money/ha]", min = 0, max = 5000, value = 0)
      )
    })
    
    
    # # Category-wise funding total for selected GOV schemes
    # category_summary <- reactive({
    #   req(input$one_schemes, input$annual_schemes)
    #   
    #   selected <- c(input$one_schemes, input$annual_schemes)
    #   
    #   filtered_data() %>%
    #     filter(`UI display` %in% selected) %>%
    #     group_by(category_name) %>%
    #     summarise(total_value = sum(Value, na.rm = TRUE)) %>%
    #     { setNames(.$total_value, .$category_name) }
    # })
    # 
    # # Overall final totals (government + private)
    # combined_totals <- reactive({
    #   cat_sum <- category_summary()
    #   
    #   private_data <- filtered_data()
    #   
    #   # Add selected one-time private funding
    #   one_private_val <- private_data %>%
    #     filter(`UI display` %in% input$one_private_schemes) %>%
    #     pull(Value) %>%
    #     sum(na.rm = TRUE)
    #   
    #   # Add selected annual private funding
    #   annual_private_val <- private_data %>%
    #     filter(`UI display` %in% input$annual_private_schemes) %>%
    #     pull(Value) %>%
    #     sum(na.rm = TRUE)
    #   
    #   cat_sum["onetime_external_support_c"] <- one_private_val + input$onetime_private_custom
    #   cat_sum["annual_external_support_c"] <- annual_private_val + input$annual_private_custom
    #   
    #   cat_sum
    #   
    #   # Display selected schemes with links
    #   output$financial_support_links <- renderUI({
    #     selected <- c(input$one_schemes,
    #                   input$annual_schemes,
    #                   input$one_private_schemes,
    #                   input$annual_private_schemes) %>% unlist() %>% unique()
    #     
    #     if (length(selected) == 0) return(NULL)
    #     
    #     selected_data <- filtered_data() %>%
    #       filter(`UI display` %in% selected) %>%
    #       select(`UI display`, Link) %>%
    #       distinct()
    #     
    #     tagList(
    #       tags$h4("Financial Support"),
    #       lapply(1:nrow(selected_data), function(i) {
    #         tags$p(
    #           tags$a(href = selected_data$Link[i], 
    #                  target = "_blank", 
    #                  rel = "noopener noreferrer",
    #                  selected_data$`UI display`[i])
    #         )
    #       })
    #     )
    #   })
    #   
    #   
    # })
    
    
    # ---------------- GOV only ----------------------------------------------
    category_summary <- reactive({
      # req(input$one_schemes || input$annual_schemes)
      
      selected_gov <- c(input$one_schemes, input$annual_schemes)
      
      cs <- 
        filtered_data() %>%                             # already country / state / Any
        filter(`UI display` %in% selected_gov) %>%
        group_by(category_name) %>%
        summarise(total_value = sum(Value, na.rm = TRUE),
                  .groups = "drop") %>%
        deframe()                                     # -> named numeric vector
      
      # print(paste0("category summary:", cs))
      cs
    })
    
    # ---------------- GOV + PRIVATE + CUSTOM --------------------------------
    combined_totals <- reactive({
      
      # start with GOV categories
      cat_sum <- category_summary()
      
      # --- add PRIVATE schemes chosen from the two dropdowns ---
      private_selected <- c(input$one_private_schemes,
                            input$annual_private_schemes) %>% unlist()
      
      if (length(private_selected) > 0) {
        priv_tbl <- filtered_data() %>%
          filter(`UI display` %in% private_selected) %>%
          group_by(category_name) %>%
          summarise(total_value = sum(Value, na.rm = TRUE),
                    .groups = "drop")
        
        # merge into cat_sum, adding if category already exists
        priv_tbl %>% purrr::pwalk(function(category_name, total_value) {
          cat_sum[category_name] <<- (cat_sum[category_name] %||% 0) + total_value
        })
      }
      
      # --- add the two custom numeric inputs into their fixed categories ----
      cat_sum["onetime_external_support_c"] <- (cat_sum["onetime_external_support_c"] %||% 0) +
        input$onetime_private_custom
      cat_sum["annual_external_support_c"]  <- (cat_sum["annual_external_support_c"]  %||% 0) +
        input$annual_private_custom
      # print(paste0("onetime_private: ", input$onetime_private_custom))
      # print(paste0("annual_private: ", input$annual_private_custom))
      # print(paste0("cat_sum: ", cat_sum))
      cat_sum          # return named numeric vector
    })
    
    # helper for %||%  (if you don't use rlang)
    `%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
    
    output$financial_support_links <- renderUI({
      selected <- c(input$one_schemes,
                    input$annual_schemes,
                    input$one_private_schemes,
                    input$annual_private_schemes) %>% 
        unlist() %>% unique()
      
      if (length(selected) == 0) return(NULL)
      
      # 1. Filter only the chosen schemes
      selected_data <- filtered_data() %>%
        filter(`UI display` %in% selected) %>%
        select(`UI display`, Link) %>%
        distinct()
      
      # 2. Split cells that contain several links (comma / semicolon / whitespace)
      library(tidyr)
      library(stringr)
      
      links_expanded <- selected_data %>%
        mutate(Link = str_split(Link, "\\s*[;,]\\s*")) %>%  # split on ; or ,
        unnest(Link)                                        # one row per link
      
      # 3. Build the UI
      tagList(
        tags$h6("For more information:"),
        lapply(split(links_expanded, links_expanded$`UI display`), function(df) {
          tags$div(
            tags$strong(df$`UI display`[1]),
            tags$ul(
              lapply(df$Link, function(url) {
                tags$li(tags$a(href = url,
                               target = "_blank",
                               rel = "noopener noreferrer",
                               url))
              })
            )
          )
        })
      )
    })
    
    return(
      list(
        country = reactive(input$country),
        state = reactive(input$state),
        selected_one_time_gov = reactive(input$one_schemes),
        selected_annual_gov = reactive(input$annual_schemes),
        onetime_private_input = reactive(input$one_private_schemes),
        annual_private_input = reactive(input$annual_private_schemes),
        category_totals = category_summary,
        total_funding_with_private = combined_totals
      )
    )
  })
}

