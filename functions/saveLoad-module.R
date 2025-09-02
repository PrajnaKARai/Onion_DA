file_dir <- if (Sys.info()[["sysname"]] == "Windows") {
  "user-states/Germany"
} else {
  "/srv/shiny-app-data/user-states/Germany"
}

# Save / Load UI ----
# ---------------------------------------------------------------------------
#  Save / Load MODULE  —  v2
# ---------------------------------------------------------------------------
saveLoadUI <- function(id) {
  ns <- NS(id)
  
  accordion_panel(
    title = "Save / Load project",
    icon  = icon("floppy-disk"),
    
    ## ---- Save area ----------------------------------------------------------
    textInput(ns("project_name"), "Project name"),
    actionButton(ns("save_btn"), "Save", icon("save"), class = "btn-primary"),
    
    ## ---- Saved versions -----------------------------------------------------
    br(),                          # just a little spacing
    selectInput(ns("version_sel"),
                "Saved versions",
                choices = character(0)),     # empty at start → will be filled
    
    div(class = "d-flex gap-2",
        actionButton(ns("load_btn"),   "Load",   icon("rotate-left")),
        actionButton(ns("delete_btn"), "Delete", icon("trash"))
    ),
    
    hr(),
    downloadButton(ns("dl_csv"), "Download current inputs (.csv)")
  )
}


# ---------------------------------------------------------------------------
#  SAVE / LOAD SERVER  –  uses current_input_table for state
# ---------------------------------------------------------------------------
saveLoadServer <- function(id, input_table_fun) {
  moduleServer(id, function(input, output, session) {
    ## ── user directory helpers ──────────────────────────────────────────────
    user_id <- reactive({
      uid <- session$user
      if (is.null(uid) || uid == "") "anon"
      else gsub("[^A-Za-z0-9_.-]", "_", uid)
    })
    
    user_dir <- reactive({
      dir <- file.path(file_dir, user_id())
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      dir
    })
    
    list_versions <- function() {
      tibble::tibble(
        file  = sort(list.files(user_dir(), pattern = "\\.rds$", full.names = TRUE)),
        label = basename(file)
      )
    }
    
    prune_to_five <- function() {
      f <- list_versions()$file
      if (length(f) > 5) unlink(head(f, -5))
    }
    
    refresh_version_choices <- function() {
      dir   <- isolate(user_dir())
      files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
      updateSelectInput(session, "version_sel",
                        choices  = basename(files),
                        selected = if (length(files)) tail(basename(files), 1)
                        else character(0))
    }
    
    observe({ refresh_version_choices() })   # fill dropdown at start
    
    ## ── grab the *current* values of *all* widgets we care about ------------
    live_state <- reactive({
     # keep every input whose id ends with _c, _p, _t, _n, or _cond
      wanted <- grep("(_c$|_p$|_t$|_n$|_cond$)", names(input), value = TRUE)
      reactiveValuesToList(input)
      reactiveValuesToList(input)[wanted]
    })
    
    ## ── smart updater (unchanged) ───────────────────────────────────────────
    update_input_smart <- function(id, value) {
      if (!id %in% names(session$input)) return()
      
      if (is.numeric(value) && length(value) == 2)
        try(updateSliderInput( session, id, value = value), silent = TRUE)
      else if (is.numeric(value) && length(value) == 1) {
        try(updateNumericInput(session, id, value = value), silent = TRUE)
        try(updateSliderInput( session, id, value = value), silent = TRUE)
      } else if (is.logical(value))
        try(updateCheckboxInput(session, id, value = value), silent = TRUE)
      else if (is.character(value) || is.factor(value)) {
        ch <- as.character(value)
        ok <- try(updateSelectInput(session, id, selected = ch), silent = TRUE)
        if (inherits(ok, "try-error"))
          try(updateTextInput(session, id, value = paste(ch, collapse = ", ")),
              silent = TRUE)
      }
    }
    
    ## ── SAVE ────────────────────────────────────────────────────────────────
    observeEvent(input$save_btn, {
      
      # print(input)
      # print(output)
      # print(session)

      
      req(nzchar(input$project_name))
      
      state <- live_state()                    # ⬅︎  real widget values
      
      fname <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S_"),
                      make.names(input$project_name), ".rds")
      saveRDS(state, file.path(user_dir(), fname))
      prune_to_five()
      refresh_version_choices()
      showNotification("Saved ✔", type = "message")
    })
    
    ## ── LOAD ────────────────────────────────────────────────────────────────
    observeEvent(input$load_btn, {
      req(input$version_sel)
      target <- file.path(user_dir(), input$version_sel)
      if (!file.exists(target)) {
        showNotification("File vanished 🤷 ♂️", type = "error"); return()
      }
      
      state <- readRDS(target)
      lapply(names(state), function(nm) update_input_smart(nm, state[[nm]]))
      showNotification("Restored ✔", type = "message")
    })
    
    ## ── DELETE ──────────────────────────────────────────────────────────────
    observeEvent(input$delete_btn, {
      req(input$version_sel)
      unlink(file.path(user_dir(), input$version_sel))
      refresh_version_choices()
      showNotification("Deleted 🗑", type = "message")
    })
    
    ## ── DOWNLOAD CSV  (unchanged) ───────────────────────────────────────────
    output$dl_csv <- downloadHandler(
      filename = function() paste0("inputs_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(input_table_fun(), file, row.names = FALSE)
    )
  })
}
