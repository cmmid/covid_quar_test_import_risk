if (!require(pacman)){
  install.packages("pacman")
}

pacman::p_load(char = c("tidyverse",
                        "here",
                        "grid",
                        "ggh4x",
                        "stringi",
                        "scales",
                        "mgcv",
                        "furrr",
                        "tictoc",
                        "conflicted",
                        "magrittr",
                        "epitools",
                        "rlang",
                        "cowplot",
                        "ggnewscale",
                        "formula.tools",
                        "measurements",
                        "fitdistrplus",
                        "magick",
                        "pdftools",
                        "truncdist",
                        "rriskDistributions",
                        "ggridges",
                        "readODS",
                        "patchwork",
                        "htmlTable",
                        "directlabels",
                        "scales",
                        "ggrepel",
                        "EnvStats",
                        "DescTools",
                        "fst",
                        "rlang",
                        "countrycode",
                        "readr","brms","broom","modelr",
                        "dtplyr",
                        "furrr",
                        "data.table",
                        "geometry",
                        "MESS",
                        "grDevices"))

conflicted::conflict_prefer("set_names", "purrr")
conflicted::conflict_prefer("melt", "reshape2")
conflicted::conflict_prefer("predict", "stats")
conflict_prefer(":=", "data.table")
map(.x = c("mutate", "select", "filter"), 
    .f = function(x){conflicted::conflict_prefer(name = x, "dplyr")})

