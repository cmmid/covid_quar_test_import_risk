# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_all.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file) %>% 
         filter(stringency=="Pre-board only") %>% 
         mutate(
           assay=fct_relevel(assay,"No testing",
                             "LFT",
                             "PCR"),
           pre_board_assay = replace_na(pre_board_assay,"No pre-flight testing"),
           pre_board_assay = fct_relevel(pre_board_assay,"No pre-flight testing",
                                         "LFT pre-flight",
                                         "PCR pre-flight")
         ) %>% 
         mutate(assay=fct_recode(assay,`LFT`="Innova (KCL)"),
                pre_board_assay=fct_recode(pre_board_assay,
                                           `LFT pre-flight`="Innova (KCL)",
                                           `PCR pre-flight`="PCR")) %>%
         mutate(adherence=ifelse(adherence_quar == 1 & 
                                   adherence_test == 1 &
                                   adherence_symp == 1, "full","literature"))%>% 
         filter(pre_board_assay!="No pre-flight testing")) 

#including symptomatic self-isolation:

pre_flight <- plotting_func(
  x             = get(results_name),
  x_var         = pre_board_test_delay,
  group_var     = pre_board_assay,
  include_symp  = F,
  my_palette    = post_test_pal,
  zero_line     = TRUE) 

pre_flight$plot <- pre_flight$plot +
  scale_color_manual(values = set_names(post_test_pal,
                                        paste(names(post_test_pal),
                                              "pre-flight")), name = "")

pre_flight_plot <- pre_flight$plot&
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())&
  labs(x="Delay from pre-flight test til flight departure (days)",
       y      = expression(Change~'in'~R[0]))&
  scale_y_continuous(labels = scales::percent)

pre_flight$summaries

save_plot(dpi    = 400,
          device = "png",
          prefix = "figS5",
          base   = "plot", 
          width  = 4.5, 
          height = 4.5,
          units  = 'inch')
