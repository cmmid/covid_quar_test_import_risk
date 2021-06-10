source("packages.R")
source("plot_functions.R")

imports_raw <- read_csv("covid_imported_cases/outputs/april_estimates.csv") %>% 
  select(-X1) %>% 
  mutate(region = countrycode(iso_code, 
                              origin = "iso3c",
                              destination="un.region.name"),
         subregion = countrycode(iso_code,
                                 origin = "iso3c",
                                 destination = "un.regionsub.name")) 

wbc <- wbstats::wbcountries('en')

# which countries are missing?
anti_join(select(filter(wbc, !is.na(capital)), 
                 iso3c, country),
          select(imports_raw, iso3c = iso_code)) %>%
  left_join(select(wbc, iso3c, region)) %>%
  mutate(region = trimws(region)) %>%
  split(.$region) %>%
  map(~select(.x, -region))


if("TWN" %in% imports_raw$iso_code) {
  imports_raw %<>% 
    rows_update(x = ., y = tibble(iso_code  = "TWN", 
                                  region    = "Asia", 
                                  subregion = "Eastern Asia"))
}


imports_raw %<>% 
  mutate(subregion = fct_recode(
    .f = subregion,
    CANZUS = "Northern America",
    CANZUS = "Australia and New Zealand",
    `Southern and Central Asia` = "Southern Asia",
    `Southern and Central Asia` = "Central Asia",
    `Latin America\nand the Caribbean`="Latin America and the Caribbean")) %>%
  mutate(region = ifelse(
    subregion == "CANZUS", "Oceania and\nNorthern America", region)) %>%
  mutate(region = factor(region, levels = c("Africa",
                                            "Americas",
                                            "Oceania and\nNorthern America",
                                            "Asia",
                                            "Europe"))) 

# imported cases under different interventions

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_redo_inf.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file))


comparison_scenarios <-
  list(`No intervention` = 
         data.frame(quar_dur        = 0,
                    pre_board_assay = "No pre-flight testing",
                    assay           = "No testing",
                    adherence_symp  = 0.71,
                    adherence_quar  = 0.28,
                    adherence_test  = 0.86),
       `Pre-flight LFT` =
         data.frame(quar_dur        = 0,
                    pre_board_test_delay = 0,
                    pre_board_assay = "LFT pre-flight",
                    assay           = "No testing",
                    adherence_symp  = 0.71,
                    adherence_quar  = 0.28,
                    adherence_test  = 0.86),
       `Pre-flight LFT + quarantine (5) + LFT` = 
         data.frame(quar_dur        = 5,
                    pre_board_test_delay = 0,
                    pre_board_assay = "LFT pre-flight",
                    assay           = "LFT",
                    adherence_symp  = 0.71,
                    adherence_quar  = 0.28,
                    adherence_test  = 0.86),
       `Pre-flight LFT + quarantine (10) + LFT` = 
         data.frame(quar_dur        = 10,
                    pre_board_test_delay = 0,
                    pre_board_assay = "LFT pre-flight",
                    assay           = "LFT",
                    adherence_symp  = 0.71,
                    adherence_quar  = 0.28,
                    adherence_test  = 0.86),
       `Pre-flight LFT + daily testing (5)` = 
         data.frame(pre_board_assay = "LFT pre-flight",
                    pre_board_test_delay = 0,
                    n_tests         = 5,
                    adherence_symp  = 0.71,
                    adherence_quar  = 0.28,
                    adherence_test  = 0.86)
  ) 


comparison_processed <- comparison_scenarios %>% 
  map(~calculate_detected_infectious_status(x = inner_join(.x, 
                                                           get(results_name)),
                                            detected_levels   = detected_levels,
                                            infectious_levels = infectious_levels)) 

# now how many were detected?
comparison_summarised <- comparison_processed %>%
  map_df(.id = "Intervention",
         ~mutate(.x,infectious_ = grepl(pattern = "^Infectious", x = infectious)) %>%
           mutate(infectious = ifelse(grepl(x = infectious,
                                            pattern = "(I|i)nfectious.*(release|adherence)"),
                                      "Risk of transmission in community",
                                      "No risk of transmission in community")) %>% 
           count(sim,
                 #detected,
                 infectious,
                 .drop = FALSE) %>%
           group_by(sim) %>%
           mutate(N = sum(n),
                  p = n/N) %>%
           ungroup %>%
           select(-n, -N) %>%
           group_by(#detected,
             infectious) %>%
           nest(data = c(sim,p)) %>%
           
           mutate(Q = map(.x = data, ~quantile(.x$p, c(0.025, 0.5, 0.975))),
                  M = map_dbl(.x = data, ~mean(.x$p))) %>%
           unnest_wider(Q) %>%
           select(-data) )

comparison_summarised_risk <- 
  filter(comparison_summarised,
         grepl(pattern = "^Risk", x = infectious))

imports <- imports_raw %>% 
  mutate(Intervention = "No intervention") %>%
  inner_join(comparison_summarised_risk) %>%
  mutate_at(.vars = vars(starts_with("imported_cases")),
            .funs = ~multiply_by(., `50%`)) %>% 
  drop_na(imported_cases_mid) %>% 
  pivot_longer(cols = starts_with('imported_cases'))


imports %<>%
  ungroup %>%
  arrange(region, subregion) %>%
  dplyr::mutate(subregion2 = factor(subregion, levels = unique(.$subregion), ordered = T)) 

value_min_y <- 0.1
value_min_x <- 1


imports_to_plot <- imports %>%
  mutate(value               = pmax(value, value_min_y)) %>%
  mutate(incidence_total_mid = pmax(incidence_total_mid, value_min_x))

set.seed(133032)
imports_plot <- imports_to_plot %>%
  filter(name == "imported_cases_mid") %>%
  ggplot(data = . )+
  geom_ribbon(data = fill_4,
              aes(x = value, ymin = ymin*value, ymax = ymax*value, fill = band),
              alpha = 0.5) +
  geom_text_repel(aes(x     = incidence_total_mid,
                      y     = value,
                      label = iso_code),
                  label.size = NA, force = 150,
                  #color = "white",
                  color = rgb(0,0,0,0.8),
                  fill  = rgb(0,0,0,0.4),
                  segment.colour = rgb(0,0,0,0.4),
                  #segment.alpha = 1,
                  #label.padding = 1,
                  #max.iter = Inf,
                  #size = 3, 
                  min.segment.length = 0,
                  max.time=20
  )+
  geom_point(aes(x    = incidence_total_mid,
                 y    = value),
             #size=2,
             shape=4)+
  scale_fill_brewer(palette = "Reds",
                    name = "Ratio of infectious arrivals to total domestic incidence") +
  coord_cartesian(ylim = c(value_min_y*0.3, 5e3),
                  expand = FALSE,
                  xlim = c(0.3 ,xrange$value[2])) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(x)label_number_si2(x,lt = TRUE)) +
  scale_x_log10(labels = function(x)label_number_si2(x,lt = TRUE),
                breaks = 10^seq(0,6,by=2)
  ) +
  labs(x = "Estimated domestic incidence per day",
       y = "Estimated number of infectious arrivals per day")+
  plotting_theme+
  theme(panel.spacing.y = unit(1, "lines")) +
  # ggpubr::rotate_x_text(angle=45)+
  # ggpubr::rotate_y_text(angle=45)+
  facet_wrap(~subregion2, ncol=4) 

imports_plot

save_plot(plot = imports_plot, prefix = "who", base = "imports_rescaled",
          device = "png", width=9, height=9, units = "inch", dpi=450)

comparisons_to_plot <- 
  comparison_summarised_risk  %>%
  arrange(desc(`50%`)) %>%
  mutate(Intervention = fct_inorder(Intervention, ordered = T)) %>%
  crossing(imports_raw %>% filter(iso_code %in% country_list)) %>%
  mutate(risk_med  = `50%` * imported_cases_mid  / incidence_total_mid,
         risk_low  = `50%` * imported_cases_low  / incidence_total_high,
         risk_high = `50%` * imported_cases_high / incidence_total_low) %>%
  select(Intervention, iso_code, country, starts_with("risk")) %>%
  mutate_at(.vars = vars(starts_with("risk")),
            .funs = ~pmax(1e-3, .)) %>%
  mutate(country = countrycode(iso_code, 'iso3c', 'country.name'),
         country2 = fct_rev(fct_reorder(country,.x = iso_code, .fun = unique)))

comparison_plot <- comparisons_to_plot %>%
  ggplot(data = ., aes(y = risk_med, x = country2)) +
  geom_ribbon(data = fill_3 %>% filter(key == "xmin",
                                       variant == FALSE) %>%
                crossing(x = c(-Inf, Inf)),
              aes(ymin = ymin*100, ymax = ymax*100,
                  alpha = 0.5,
                  x = x,
                  fill = band), show.legend = FALSE,
              inherit.aes = FALSE) +
  scale_fill_brewer(palette="Reds")+
  geom_linerange(aes(ymin = risk_low,
                     ymax = risk_high,
                     group = rev(Intervention)),
                 position = position_dodge(width = 0.75)) +
  geom_point(aes(shape = Intervention,
                 group = rev(Intervention)
  ),
  fill="white",
  size=2,
  position = position_dodge(width = 0.75)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = axis_label_lt_percent,
                limits = c(NA, 11)) + 
  coord_flip() +
  plotting_theme +
  labs(x="", y="Ratio of undetected/non-adherent infectious arrivals to total domestic incidence") + 
  guides(shape=guide_legend(ncol=3, byrow=F,reverse=F)) +
  #theme(legend.position = 'right') +
  scale_shape_manual(values = shapes, breaks = rev(names(shapes)), name = NULL) 


save_plot(plot = comparison_plot, prefix = "who", base = "intervention_comparison",
          device = "png", width=9, height=4.5, units = "inch", dpi=450)


comparisons_to_plot_all_countries <- 
  comparison_summarised_risk  %>%
  arrange(desc(`50%`)) %>%
  mutate(Intervention = fct_inorder(Intervention, ordered = T)) %>%
  crossing(imports_raw %>% arrange(region, subregion) %>%
             dplyr::mutate(subregion2 = factor(subregion, 
                                               levels = unique(.$subregion),
                                               ordered = T))) %>%
  mutate(risk_med  = `50%` * imported_cases_mid  / incidence_total_mid,
         risk_low  = `50%` * imported_cases_low  / incidence_total_high,
         risk_high = `50%` * imported_cases_high / incidence_total_low) %>%
  select(Intervention, iso_code, country, subregion2, starts_with("risk")) %>%
  mutate_at(.vars = vars(starts_with("risk")),
            .funs = ~pmax(1e-4, .)) %>%
  mutate(country = countrycode(iso_code, 'iso3c', 'country.name'),
         country2 = fct_rev(fct_reorder(country,.x = iso_code, .fun = unique))) %>% 
  drop_na(risk_med)

plot_all_countries <- comparisons_to_plot_all_countries %>%
  nest(data = -c(subregion2)) %>% 
  mutate(id=row_number()) %>% 
  mutate(plot_=map2(data,subregion2,~ggplot(data = ., aes(y = risk_med, x = country2)) +
                      geom_ribbon(data = fill_3 %>% filter(key == "xmin",
                                                           variant == FALSE) %>%
                                    crossing(x = c(-Inf, Inf)),
                                  aes(ymin = ymin, ymax = ymax,
                                      alpha = 0.5,
                                      x = x,
                                      fill = band), show.legend = FALSE,
                                  inherit.aes = FALSE) +
                      scale_fill_brewer(palette="Reds")+
                      geom_linerange(aes(ymin = risk_low,
                                         ymax = risk_high,
                                         group = rev(Intervention)),
                                     position = position_dodge(width = 0.75)) +
                      geom_point(aes(shape = Intervention,
                                     group = rev(Intervention)
                      ),
                      fill="white",
                      size=2,
                      position = position_dodge(width = 0.75)) +
                      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = axis_label_lt_percent,
                                    limits = c(NA, 11)) + 
                      coord_flip() +
                      plotting_theme +
                      labs(x="", y="Ratio of undetected/non-adherent infectious arrivals to total domestic incidence") + 
                      guides(shape=guide_legend(ncol=3, byrow=F,reverse=F)) +
                      #theme(legend.position = 'right') +
                      scale_shape_manual(values = shapes, breaks = rev(names(shapes)), name = NULL)+
                      ggtitle(.y)
  ))

align_all_plots <- patchwork::align_plots(c(plot_all_countries$plot_))

purrr::safely(pmap(list(plot = align_all_plots, 
                        prefix = "who", 
                        base = paste0("intervention_comparison_all_",
                                      plot_all_countries$id),
                        device = "png", 
                        width=9, 
                        height=4.5, 
                        units = "inch", 
                        dpi=450),
                   .f=save_plot))

