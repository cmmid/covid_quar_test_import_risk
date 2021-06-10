# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_inf.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file) %>% 
         filter(stringency!="Pre-board only") %>% 
         mutate(
           assay=fct_relevel(assay,"No testing",
                             "LFT",
                             "PCR"),
           pre_board_assay = replace_na(pre_board_assay,"No pre-flight testing"),
           pre_board_assay = fct_relevel(pre_board_assay,"No pre-flight testing",
                                         "LFT pre-flight",
                                         "PCR pre-flight")
         ) %>% 
         # mutate(assay=fct_recode(assay,`LFT`="Innova (KCL)"),
         #        pre_board_assay=fct_recode(pre_board_assay,
         #                                   `LFT pre-flight`="Innova (KCL)",
         #                                   `PCR pre-flight`="PCR")) %>% 
         mutate(adherence=ifelse(adherence_quar == 1 & 
                                   adherence_test == 1 &
                                   adherence_symp == 1, "full","literature"))) 

#including symptomatic self-isolation:

inc_symp_a <- plotting_func(x             = get(results_name),
                            x_var         = quar_dur,
                            faceting      = pre_board_assay ~ adherence,
                            group_var     = assay,
                            include_symp  = T, 
                            my_palette    = post_test_pal,
                            zero_line     = FALSE)

inc_symp_b <- plotting_func(x             = get(results_name),
                            x_var         = n_tests,
                            faceting      = pre_board_assay ~ adherence,
                            group_var     = assay,
                            include_symp  = T, 
                            my_palette    = post_test_pal,
                            zero_line     = FALSE)

inc_symp_a$summaries %>% filter(adherence=="literature") %>% ungroup()%>% slice(1)

figS3 <- inc_symp_a$plot + (inc_symp_b$plot+guides(colour=F)) +
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(5,4), guides = "collect")&
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) &
  labs(y      = expression(Change~'in'~R[0]),
       colour = "Post-flight test") &
  scale_y_continuous(labels = scales::percent, limits = c(NA, -0.25))

save_plot(dpi    = 400,
          plot   = figS3,
          device = "png",
          prefix = "figS3",
          base   = "plot", 
          width  = 9, 
          height = 4.5,
          units  = 'inch')


a <- plotting_func(x             = get(results_name),
                   x_var         = quar_dur,
                   faceting      = pre_board_assay ~ adherence,
                   group_var     = assay,
                   my_palette    = post_test_pal,
                   zero_line     = TRUE)

a$summaries %>% filter(adherence=="literature") 
a$summaries %>% filter(adherence=="literature") %>% filter(quar_dur==14)
a$summaries %>% filter(adherence=="literature") %>% filter(quar_dur==5)

a$summaries %>% filter(adherence=="full") %>% filter(quar_dur%in%c(10,14)) %>%
  split(.$pre_board_assay == "No pre-flight testing")

b <- plotting_func(x             = get(results_name),
                   x_var         = n_tests,
                   faceting      = pre_board_assay ~ adherence,
                   group_var     = assay,
                   my_palette    = post_test_pal,
                   zero_line     = TRUE)


b$summaries %>% filter(adherence=="literature") %>% 
  filter(n_tests%in%c(3,5,7,10)) %>%
  split(.$pre_board_assay)

fig6 <- a$plot+(b$plot+guides(colour=F))+
  plot_annotation(tag_levels = "A") + 
  plot_layout(widths = c(5,4),guides = "collect") &
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) &
  labs(y= expression(Change~'in'~R[s]),
       colour="Post-flight test") &
  scale_y_continuous(labels = scales::percent,limits = c(NA,0.1))

save_plot(dpi    = 400,
          plot   = fig6,
          device = "png",
          prefix = "fig6",
          base   = "plot", 
          width  = 9, 
          height = 4.5,
          units  = 'inch')