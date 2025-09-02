txt <- readLines("functions/Walnut_grain_veg_tub_mcsim-only.R")
# prepend every symbol that ends in _c/_p/_t/_n with x$
txt <- gsub("\\b([A-Za-z0-9_]+_[cptn])\\b", "x[['\\1']]", txt, perl = TRUE)
writeLines(txt, "functions/Walnut_grain_veg_tub_df.R")
