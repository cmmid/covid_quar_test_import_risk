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
  labs(y= expression(Change~'in'~R),
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

#####################
#### RISK RATIOS ####
#####################

baseline_low <- get(results_name) %>%
  as_tibble() %>% 
  filter(quar_dur  == 0,
         pre_board_assay==c("No pre-flight testing"),
         adherence_quar!=1,
         adherence_test!=1,
         adherence_symp!=1,
         assay=="No testing") %>% 
  select(sim, max_auc,tot_auc, flight_arr_t) %>% 
  group_by(sim) %>% 
  summarise(baseline_prop=1-(sum(max_auc)/sum(tot_auc))) %>% 
  mutate(name="baseline_low")

baseline_high <- get(results_name) %>%
  as_tibble() %>% 
  filter(quar_dur  == 14,
         pre_board_assay==c("No pre-flight testing"),
         adherence_quar!=1,
         adherence_test!=1,
         adherence_symp!=1,
         assay=="No testing") %>% 
  select(sim, max_auc,tot_auc, flight_arr_t) %>% 
  group_by(sim) %>% 
  summarise(baseline_prop=1-(sum(max_auc)/sum(tot_auc))) %>% 
  mutate(name="baseline_high")

#### vs low baseline ----
a <- get(results_name) %>% 
  filter(
    pre_board_assay==c("No pre-flight testing"),
    adherence_quar!=1,
    adherence_test!=1,
    adherence_symp!=1) %>% 
  # needs to be harmonious with the plotting_func()
  rr_func(x=.,
          baseline = baseline_low,
          x_var=quar_dur,
          group_var=assay,
          log = T)

b <- get(results_name) %>% 
  filter(
    pre_board_assay==c("No pre-flight testing"),
    adherence_quar!=1,
    adherence_test!=1,
    adherence_symp!=1) %>% 
  rr_func(x=.,
          baseline = baseline_low,
          x_var=n_tests,
          group_var  =assay,
          log = T) 

fig3 <- a$plot+(b$plot+guides(colour=FALSE)+labs(y=""))

p_ranges_y <- c(10^(ggplot_build(fig3[[1]])$layout$panel_scales_y[[1]]$range$range),
                10^(ggplot_build(fig3[[2]])$layout$panel_scales_y[[1]]$range$range))

fig3+plot_annotation(tag_levels = "A")+
  plot_layout(widths = c(3,2),guides = "collect")&
  theme(legend.position = "bottom")&
  scale_y_continuous(name = "Relative change in R compared to baseline",limits = c(min(p_ranges_y), max(p_ranges_y)))

save_plot(dpi = 400, 
          device = "png",
          prefix = "baseline_low_R",
          base = "RR_plot", 
          width = 210, 
          height = 120)

#### adherence -----
a <- get(results_name) %>% 
  filter(
    pre_board_assay==c("No pre-flight testing")) %>% 
  rr_func(x=.,
          baseline = baseline_low,
          x_var=quar_dur,
          row_vars = adhering_iso,
          col_vars = adhering_quar,
          group_var=assay,
          log = T)

b <- get(results_name) %>% 
  filter(
    pre_board_assay==c("No pre-flight testing")) %>% 
  rr_func(x=.,
          baseline = baseline_low,
          x_var=n_tests,
          row_vars = adhering_iso,
          group_var  =assay,
          log = T) 

fig4 <- a$plot+(b$plot+guides(colour=FALSE)+labs(y=""))

p_ranges_y <- c(10^(ggplot_build(fig4[[1]])$layout$panel_scales_y[[1]]$range$range),
                10^(ggplot_build(fig4[[2]])$layout$panel_scales_y[[1]]$range$range))

fig4+plot_annotation(tag_levels = "A")+
  plot_layout(widths = c(3,2),guides = "collect")&
  theme(legend.position = "bottom")&
  scale_y_continuous(name = "Relative change in R compared to baseline",limits = c(min(p_ranges_y), max(p_ranges_y)))

save_plot(dpi = 400, 
          device = "png",
          prefix = "baseline_low_adherence",
          base = "RR_plot", 
          width = 210, 
          height = 120)


##### vs. 14 day quarantine ----
a <- get(results_name) %>% 
  #filter(pre_board_assay==c("No pre-flight testing")) %>% 
  rr_func(x=.,
          baseline = baseline_high,
          x_var=quar_dur,
          col_vars = pre_board_assay,
          row_vars = adherence,
          group_var=assay,
          log = T)

b <- get(results_name) %>% 
  #filter(pre_board_assay==c("No pre-flight testing")) %>% 
  rr_func(x=.,
          baseline = baseline_high,
          x_var=n_tests,
          col_vars = pre_board_assay,
          row_vars = adherence,
          group_var  =assay,
          log = T) 

fig4 <- a$plot+(b$plot+guides(colour=FALSE)+labs(y=""))

p_ranges_y <- c(10^(ggplot_build(fig4[[1]])$layout$panel_scales_y[[1]]$range$range),
                10^(ggplot_build(fig4[[2]])$layout$panel_scales_y[[1]]$range$range))

fig4+plot_annotation(tag_levels = "A")+
  plot_layout(widths = c(3,2),guides = "collect")&
  theme(legend.position = "bottom")&
  scale_y_continuous(name = "Relative change in R compared to baseline",limits = c(min(p_ranges_y), max(p_ranges_y)))

save_plot(dpi = 400, 
          device = "png",
          prefix = "baseline_low_adherence",
          base = "RR_plot", 
          width = 210, 
          height = 120)