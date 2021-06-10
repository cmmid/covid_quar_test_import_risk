library(pacman)
pacman::p_load(char = c("data.table",
                        "forcats",
                        "tidyverse",
                        "ggrepel",
                        "scales",
                        "wbstats",
                        "latex2exp",
                        "stringr",
                        "countrycode",
                        #"eurostat",
                        "wbstats"))

source("plot_functions.R")


# country_list <- eu_countries %>%
#   dplyr::filter(code != "UK") %>%
#   mutate(iso_a3 = countrycode::countrycode(name, origin = "country.name", destination = "iso3c")) %>%
#   pull(iso_a3)

month.labs <- c("October", "November", "December", "January", "February", "March", "April") %>% 
  rlang::set_names(., str_to_lower(.))

month.labs %>%
  set_names(names(.),.) %>%
  as.list %>%
  map(~paste0("covid_imported_cases/outputs/", 
              paste0("uk_exports_", paste0(.x), ".csv"))) %>%
  map_df(fread, .id = "month") %>%
  distinct(destination_country_iso_code, month, incidence_mid)

exports <- fread("covid_imported_cases/outputs/uk_exports_all.csv") %>%
  mutate(year = ifelse(month %in% c("october", "november", "december"),
                       2020,
                       2021))

wbc <- wbcountries('en')

set.seed(2020)
# exports %>%
#   na.omit %>%
#   distinct(iso_a3, month) %>%
#   count(iso_a3) %>%
#   filter(n >= 6) %>%
#   mutate(iso3c = iso_a3) %>%
#   filter(!(iso3c %in% c("GBR", "IND"))) %>%
#   left_join(select(wbc, iso3c, incomeID)) %>%
#   filter(incomeID %in% c("HIC", "UMC")) %>%
#   na.omit %>%
#   group_by(incomeID) %>%
#   sample_n(2) %>% pull(iso3c) %>% c("USA", "SGP") -> country_list
# 
 
# #select(-X1) %>% 
# mutate(region = countrycode(iso_code, 
#                             origin = "iso3c",
#                             destination="un.region.name"),
#        subregion = countrycode(iso_code,
#                                origin = "iso3c",
#                                destination = "un.regionsub.name"))

# need to do better region names
# 
# source("find_latest_results.R")
# assign(results_name, read.fst(here::here("results", results_file)))


exports_to_plot <- exports %>% 
  as_tibble %>%
  filter(iso_a3 %in% country_list) %>%
  crossing(data.frame(Intervention = c("No intervention",
                                       "Pre-flight LFT + quarantine (5) + LFT"))) %>%
  inner_join(comparison_summarised_risk) %>%
  mutate_at(.vars = vars(starts_with("exports")),
            .funs = ~multiply_by(., `50%`)) %>% 
  drop_na(exports_mid) %>% 
  #pivot_longer(cols = starts_with('exports')) %>%
  mutate(month = factor(month, 
                        levels = names(month.labs),
                        labels = month.labs)) 


variant.labs <- c("Non-B.1.1.7", "B.1.1.7")
names(variant.labs) <- c(FALSE, TRUE)

exports_to_plot %<>% as.data.table %>%
  arrange(year, month) %>%
  mutate(date = paste(month, year)) %>%
  mutate(date = fct_inorder(factor(date), ordered = T)) %>%
  mutate_at(.vars = vars(contains("exports")),
            .funs = ~pmax(0.1, .)) %>%
  mutate_at(.vars = vars(starts_with("incidence")),
            .funs = ~pmax(0.1, .)) %>%
  filter(iso_a3 %in% country_list) %>%
  as.data.table 

exports_plot_time <- 
  exports_to_plot %>%
  as_tibble %>%
  dplyr::select(-contains("low"), -contains("high")) %>%
  #gather(key, value, intervention_exports_mid, exports_mid) %>%
  mutate(Risk = exports_mid/incidence_mid) %>%
  dplyr::select(-contains("mid")) %>%
  mutate(variant = variant.labs[variant + 1]) %>%
  mutate(Risk = pmax(Risk, 1e-1 * 1e-2)) %>%
  mutate(date = paste(substr(month, 1, 3), year),
         date = fct_inorder(factor(date), ordered = T)) %>%
  ggplot(data = .,
         aes(x = date, y = Risk)) +
  geom_ribbon(data = fill_4 %>% 
                crossing(x = c(-Inf, Inf)),
              aes(ymin = ymin, ymax = ymax,
                  alpha = 0.5,
                  x = x,
              fill = band), show.legend = FALSE,
              inherit.aes = FALSE) +
  scale_fill_brewer(palette="Blues")+
  ggnewscale::new_scale_fill()+
  geom_path(aes(group = interaction(date,variant)),
            position = position_dodge(width = 0.5),
            linetype="solid")+
  geom_point(aes(group = variant,
                 shape = Intervention,
                 fill = variant),
             size=2,
             position = position_dodge(width = 0.5)) +
  facet_wrap(~ iso_a3, nrow = 2,
             labeller = labeller(iso_a3 = function(x){
               countrycode::countrycode(x, "iso3c", "country.name")})) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = axis_label_lt_percent) +
  plotting_theme +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  scale_fill_manual(name = "Variant", values = c("black","white"),
                    breaks=rev(variant.labs)) +
  scale_shape_manual(name = "Scenario", 
                     values = shapes, breaks = names(shapes)) +
    guides(fill = guide_legend(override.aes = list(fill = c("black","white"),
                                                   shape=c(21,21))))+
  #annotation_logticks(sides = 'lr') +
  labs(x = "Date",
       y = "Undetected/non-adherent infectious arrivals\nfrom the UK as a percentage of total domestic incidence")


save_plot(plot = exports_plot_time, prefix = "who", base = "exports_plot_time_",
          device = "png", width=9, height=6, units = "inch", dpi=450)


# how big is each country's epidemic over time?

mutate(as.data.frame(exports_to_plot),
       risk_mid  = exports_mid/incidence_mid,
       risk_low  = exports_low/incidence_high,
       risk_high = exports_high/incidence_low) %>%
  mutate(risk = sprintf("%0.0f%% (%0.0f%%, %0.0f%%)", 
                        100*risk_mid, 
                        100*risk_low, 
                        100*risk_high)) %>%
  select(variant, date, risk, iso_a3, Intervention) %>%
  inner_join(exports_to_plot %>%
               mutate(incidence = sprintf("%s (%s, %s)", 
                                          round(incidence_mid),
                                          round(incidence_low),
                                          round(incidence_high))) %>%
               select(date, incidence, iso_a3), copy = T) %>%
  distinct %>%
  spread(variant, risk) %>%
  rename(`Risk: Non-B.1.1.7` = `FALSE`,
         `Risk: B-1.1.7`     = `TRUE`) %>%
  arrange(iso_a3, date) %>%
  mutate(Country = countrycode::countrycode(iso_a3, "iso3c", "country.name.en")) %>%
  select(Country, Month = date, Incidence = incidence, everything(), -iso_a3) %>%
  write_csv("results/incidence_exportation.csv")


