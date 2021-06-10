
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")
source("lfa_test_data.R")

run_model <- function(
  input,
  trajectories,
  seed             = 145,
  asymp_parms
){
  
  set.seed(seed)
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))

# generate testing times
my_message("Calculating test times")
#browser()

trajectories %<>% 
  crossing(distinct(input)) %>%     
  mutate(flight_arr_t  = flight_dep_t + dur_flight) %>% 
  mutate(test_t = pmap(.l = list(flight_dep_t=flight_dep_t,
                                 flight_arr_t=flight_arr_t,
                                 multiple_tests=multiple_tests,
                                 tests=tests,
                                 sampling_freq=sampling_freq,
                                 quar_dur=quar_dur,
                                 n_tests=n_tests),
                      .f = test_times)) %>% 
  unnest(test_t) 


#calc outcomes 
my_message("Calculating outcomes for each traveller")
trajectories %<>% calc_outcomes() 

#find earliest positive test result
#browser()
trajectories %<>% 
  nest(test_t, test_p, test_no, test_label, ct, screen) %>% 
  mutate(earliest_t      =  map(.f = earliest_pos, 
                                .x = data)) %>% 
  unnest_wider(earliest_t) %>% 
  rename("earliest_t"=test_t)%>% 
  select(-data)

#pre flight screening (tests and symptoms)
trajectories %<>% calc_outcomes_pre()

#shift other timings relative to onset
trajectories %<>%
  mutate(quar_end_t     = flight_arr_t + quar_dur,
         symp_end_t     = onset_t      + post_symptom_window,
         test_iso_end_t = earliest_t   + post_symptom_window)

# calculate remaining transmission potential averted by positive test

my_message("Calculating remaining transmission potential for each traveller")
#browser()
trajectories %<>% calc_auc()
 
return(trajectories)

}

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list( 
    `Daily testing` =
      crossing(sampling_freq    = 1,
               tests            = TRUE,
               multiple_tests   = TRUE,
               n_tests          = c(3,5,7,10),
               assay            = c("Innova (KCL)"),
               quar_dur         = NA,
               adherence_quar      = c(0.28,1),
               adherence_symp      = c(0.71,1),
               adherence_test      = c(0.86,1),
               pre_board_assay     = c(NA,"Innova (KCL)","PCR"),
               pre_board_test_delay = 0),
    `Pre-board only` =
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA,
               quar_dur         = NA,
               adherence_quar      = c(1),
               adherence_symp      = c(1),
               adherence_test      = c(1),
               pre_board_test_delay = c(0,3,5,7,10),
               pre_board_assay     = c(NA,"Innova (KCL)","PCR")),
    `Post-flight quarantine only` =
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA,
               quar_dur         = c(0, 5, 7, 10, 14),
               adherence_quar      = c(0.28,1),
               adherence_symp      = c(0.71,1),
               adherence_test      = c(0.86,1),
               pre_board_test_delay = c(0),
               pre_board_assay     = c(NA,"Innova (KCL)","PCR")),
    `Post-flight quarantine with LFA test` =
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               assay            = "Innova (KCL)",
               quar_dur         = c(0, 5, 7, 10, 14),
               adherence_quar      = c(0.28,1),
               adherence_symp      = c(0.71,1),
               adherence_test      = c(0.86,1),
               pre_board_assay     = c(NA,"Innova (KCL)","PCR"),
               pre_board_test_delay = 0),
    `Post-flight quarantine with PCR test` =
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               assay            = "PCR",
               quar_dur         = c(0, 5, 7, 10, 14),
               adherence_quar      = c(0.28,1),
               adherence_symp      = c(0.71,1),
               adherence_test      = c(0.86,1),
               pre_board_assay     = c(NA,"Innova (KCL)","PCR"),
               pre_board_test_delay = 0)
    ) %>%
      bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window = 10,
           synd_screen_prob    = 0,
           dur_flight          = 6/24) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>%
  rowwise %>%
  group_split()

#Create individuals with viral load trajectories
trajectories <- make_trajectories(n_cases     = 1000,
                                  n_sims      = 1000,
                                  seed        = 1000,
                                  input       = input,
                                  asymp_parms = asymp_fraction)

results_name <- "results_df"

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           trajectories=trajectories,
           seed = 1000
         )))

st=format(Sys.time(), "%Y%m%d_%H%M%S")
write.fst(get(results_name) %>% 
            map(as_tibble) %>% 
            bind_rows() %>% 
            select(-c(rx,ry,m,inf_curve)) %>% 
            as.data.frame(),paste0("results/results_",st,"_all.fst"))
saveRDS(trajectories,paste0("results/traj_",st,"_all.rds"))


