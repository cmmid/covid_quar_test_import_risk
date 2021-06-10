source("packages.R")
source("utils.R")
source("plot_functions.R")
source("find_latest_results.R")

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_redo_inf.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file) %>% 
         filter(stringency!="Pre-board only") 
       # mutate(
       #   assay = factor(assay, levels = c("Innova (KCL)", "PCR")),
       #   assay = fct_explicit_na(assay, "No testing"),
       #   pre_board_assay=replace_na(pre_board_assay,"No pre-flight testing"),
       #   pre_board_assay=fct_relevel(pre_board_assay,"No pre-flight testing",
       #                               "Innova (KCL)",
       #                               "PCR")
       # ) %>% 
       # mutate(assay = fct_recode(assay,`LFT`="Innova (KCL)"),
       #        pre_board_assay=fct_recode(pre_board_assay,
       #                                   `LFT pre-flight`="Innova (KCL)",
       #                                   `PCR pre-flight`="PCR")) %>% 
       # mutate(adherence = ifelse(adherence_quar   == 1 & 
       #                             adherence_test == 1 & 
       #                             adherence_symp == 1,
       #                           "full", "literature"))
) 

processed <- calculate_detected_infectious_status(x                 = results_df,
                                                  detected_levels   = detected_levels,
                                                  infectious_levels = infectious_levels) 

# appendix version of plot
# takes forever to run
# make_results_plot(
#   processed,
#   faceting = assay + adherence_quar + adherence_symp + adherence_test ~ pre_board_assay,
#                   prefix = "who", base = "appendix",
#                   width=9, height=18, device = "png", units = "inch")

# focus on main
baseline_scenarios <- crossing(data.frame(adherence_symp  = c(0.71, 1),#0.67,
                                          adherence_quar  = c(0.28, 1),#0.5,
                                          adherence_test  = c(0.86, 1)),
                               assay           = c("No testing", "PCR", "LFT"),
                               pre_board_assay = c("No pre-flight testing",
                                                   "PCR pre-flight", 
                                                   "LFT pre-flight")) %>%
  mutate(assay           = factor(assay, 
                                  #ordered = T,
                                  levels = c("No testing", "LFT", "PCR")),
         pre_board_assay = factor(pre_board_assay, 
                                  #ordered = T,
                                  levels = c("No pre-flight testing", 
                                             "LFT pre-flight", 
                                             "PCR pre-flight")))

processed_df <- mutate(processed,
                       assay           = factor(assay, levels = levels(assay), ordered = F),
                       pre_board_assay = factor(pre_board_assay, 
                                                ordered = F,
                                                levels = c("No pre-flight testing", 
                                                           "LFT pre-flight", 
                                                           "PCR pre-flight"))) %>%
  inner_join(baseline_scenarios) %>% 
  ungroup() %>% 
  filter(detected %!in% c("Ceased being infectious before flight",
                          "Develops symptoms before flight")) %>% 
  droplevels()

make_results_plot(x        = processed_df,
                  faceting = pre_board_assay + adherence_quar + adherence_symp + 
                    adherence_test ~ assay,
                  prefix   = "who",
                  base     = "daily_quar_baseline_and_ideal",
                  width    = 9,
                  height   = 12, 
                  device   = "png", 
                  units    = "inch")

make_results_plot(x        = processed_df %>%
                    filter(adherence_symp == 0.71,
                           adherence_quar == 0.28,
                           adherence_test == 0.86),
                  faceting = pre_board_assay  ~ assay,
                  prefix   = "who",
                  base     = "daily_quar_baseline",
                  width    = 9,
                  height   = 7, 
                  device   = "png", 
                  units    = "inch")


# focus on main but with both LFA sensitivities
# is this no longer needed?
# make_results_plot(x = processed %>% inner_join(dplyr::select(baseline_scenario, -sens_LFA)),
#                   faceting = sens_LFA ~ stringency + pre_board_test,
#                   prefix = "who", base = "daily_quar_both_sens",
#                   width=10, height=5, device = "png", units = "inch")



make_results_plot(x        = processed_df %>%
                    filter(adherence_symp == 1,
                           adherence_quar == 1,
                           adherence_test == 1),
                  faceting = pre_board_assay + type ~ assay,
                  prefix   = "who", 
                  base     = "daily_quar_both_type",
                  width    = 9, 
                  height   = 9, 
                  device   = "png",
                  units    = "inch")

# make_results_plot(x = processed %>% inner_join(dplyr::select(baseline_scenario, -sens_LFA)),
#                   faceting = sens_LFA + type ~ stringency + pre_board_test,
#                   prefix = "who", base = "daily_quar_both_sens_type",
#                   width=10, height=6, device = "png", units = "inch")
