# ============================================================================ #
#  dynamic-helper.R – build UI elements from Walnut_grain_veg_tub.xlsx rows    #
#  v2: supports header4 / horizontal line / break rows                         #
# ============================================================================ #

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
sanitize_id <- function(x) gsub("[^A-Za-z0-9]", "_", x)

create_ui_element <- function(row) {
  
  # ---------- quick check for layout-only rows -------------------------------
  special <- tolower(trimws(row[["ui_type"]]))
  if (special %in% c("header4", "horizontal line", "break")) {
    
    ui_obj <- switch(
      special,
      "header4"        = h4(row[["name"]] %||% ""),
      "horizontal line"= tags$hr(),
      "break"          = tags$br()
    )
    
    # still wrap in category filter if column exists --------------------------
    category <- trimws(row[["Expertise"]] %||% "")
    if (category == "" || is.na(category)) return(ui_obj)
    
    cats_vec <- trimws(unlist(strsplit(category, ";|,")))
    if (length(cats_vec) == 0) return(ui_obj)
    
    cat_inputs <- sprintf("input['cat_%s']", sanitize_id(cats_vec))
    show_all   <- paste0(
      "Object.keys(input).filter(k=>k.startsWith('cat_')).",
      "every(k=>input[k]===false)"
    )
    js_cond <- sprintf("(%s) || (%s)", show_all,
                       paste(cat_inputs, collapse = " || "))
    return(conditionalPanel(js_cond, ui_obj))
  }
  
  # ---------- MODIFIED logic for real input widgets --------------------------
  input_id      <- row[["variable"]]
  input_alias   <- row[["name"]]
  label         <- row[["description"]]
  min_val       <- as.numeric(row[["lower"]])
  max_val       <- as.numeric(row[["upper"]])
  distr         <- row[["distribution"]]
  ui_type       <- row[["ui_type"]]
  ui_step       <- as.numeric(row[["ui_steps"]])
  ui_cond       <- row[["ui_conditional"]]
  ui_cond_nam   <- row[["ui_conditional_name"]]
  ui_opt_nam    <- row[["ui_option_name"]]
  ui_opt_var    <- row[["ui_option_variable"]]
  category      <- row[["Expertise"]] %||% ""
  crop_rotation <- row[["Crop_rotation"]] %||% ""
  
  default        <- min_val
  default_2side  <- c(min_val, max_val)
  if  (distr == "tnorm_0_1") default_2side  <- c(0.05, 0.95)
  
  if (distr == "posnorm") {
    # min_val <- 0
    min_val <- min_val - abs(min_val)*0.5
    if (min_val<=0) min_val <- 0.0001
  } else if  (distr == "tnorm_0_1"){
    min_val <- 0
  } else min_val <- min_val - abs(min_val)*0.5
  
  if (distr == "tnorm_0_1") {
    max_val <- 1
  } else max_val <- max_val + abs(max_val)*0.5
  
  ui_cond <- if (!is.null(ui_cond)) as.logical(ui_cond) else FALSE
  
  if (startsWith(ui_type, "select")) {
    ui_opt_nam  <- unlist(strsplit(ui_opt_nam,  " next_option "))
    ui_opt_var  <- unlist(strsplit(ui_opt_var,  " next_option "))
    names(ui_opt_var) <- ui_opt_nam
  }
  
  # tooltip_wrapper <- function(ui_element, tooltip_text) {
  #   withTags(div(title = tooltip_text, ui_element))
  # }
  
  
  label_with_tooltip <- function(input_id, input_alias, tooltip) {
    tags$label(`for` = input_id,
               title = tooltip,
               style = "font-weight: bold; cursor: help;",
               input_alias)
  }
  
  real_ui <- switch(
    ui_type,
    "header"  = tagList(
      tags$span(style="font-size:120%;font-weight:bold;", label),
      br(), br()),
    "slider1" = tagList(
      label_with_tooltip(input_id, input_alias, label),
      sliderInput(input_id, NULL, min_val, max_val, default, step = ui_step)
    ),
    "slider2" = tagList(
      label_with_tooltip(input_id, input_alias, label),
      sliderInput(input_id, NULL, min_val, max_val, default_2side, step = ui_step)
    ),
    "numeric" = tagList(
      label_with_tooltip(input_id, input_alias, label),
      numericInput(input_id, NULL, default, step = ui_step)
    ),
    "select"  = tagList(
      label_with_tooltip(input_id, input_alias, label),
      selectInput(input_id, NULL, choices = ui_opt_var)
    ),
    "select2" = tagList(
      label_with_tooltip(input_id, input_alias, label),
      selectInput(input_id, NULL, choices = ui_opt_var, multiple = TRUE)
    ),
    "checkbox" = tagList(
      checkboxInput(
        input_id,
        label_with_tooltip(input_id, input_alias, label),
        value = as.logical(default)   # default ← lower column or an explicit “default” field
   )
   ),
    # fallback
    tagList(
      label_with_tooltip(input_id, input_alias, label),
      textInput(input_id, NULL, "")
    )
  )
  
  if (ui_cond) {
    toggle_id <- paste0(input_id, "_toggle")
    real_ui <- tagList(
      checkboxInput(toggle_id, ui_cond_nam, FALSE),
      conditionalPanel(sprintf("input['%s']", toggle_id), real_ui)
    )
  }
  
  # ---------- category filter wrapper for normal elements --------------------
  category <- trimws(category)
  crop_rotation <- trimws(crop_rotation)
  
  if ((category == "" || is.na(category)) &&
      (crop_rotation == "" || is.na(crop_rotation))) return(real_ui)
  
  cats_vec <- trimws(unlist(strsplit(category, ";|,")))
  crop_vec <- trimws(unlist(strsplit(crop_rotation, ";|,")))
  crop_vec <- crop_vec[crop_vec != "" & !is.na(crop_vec) &
                         tolower(crop_vec) != "na"]
  
  if (length(cats_vec) == 0 && length(crop_vec) == 0) return(real_ui)
  
  # ---------------------------------------------------------------------------
  # A.  Expertise logic  (opt-out: show rows when no cat boxes ticked)
  # ---------------------------------------------------------------------------
  if (length(cats_vec)) {
    cat_inputs <- sprintf("input['cat_%s']",
                          sanitize_id(cats_vec))
    
    cat_show_all <- paste0(
      "Object.keys(input).filter(k => k.startsWith('cat_')).",
      "every(k => input[k] === false)"
    )
    
    cat_cond <- sprintf("(%s) || (%s)",
                        cat_show_all,
                        paste(cat_inputs, collapse = ' || '))
  } else {
    cat_cond <- "true"          # row has no Expertise → ignore that filter
  }
  
  # ---------------------------------------------------------------------------
  # B.  Crop-rotation logic  (opt-in: hide rows with crops until selected)
  # ---------------------------------------------------------------------------
  if (length(crop_vec)) {
    crop_inputs <- sprintf("input['crop_%s']",
                           sanitize_id(crop_vec))
    crop_cond <- paste(crop_inputs, collapse = ' || ')
  } else {
    crop_cond <- "true"         # row has no crop tag → always pass
  }
  
  # ---------------------------------------------------------------------------
  # C.  Combine the two
  # ---------------------------------------------------------------------------
  js_cond <- sprintf("(%s) && (%s)", cat_cond, crop_cond)
  
  conditionalPanel(js_cond, real_ui)
  
}
