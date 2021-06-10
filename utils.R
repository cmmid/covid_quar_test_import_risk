my_message <- function(x, ...){
  message(paste(Sys.time(), x, sep = "    "), ...)
}

plan(multisession,workers=6)
options(future.globals.maxSize= 2000*1024^2)

probs        <- c(0.025,0.25,0.5,0.75,0.975)

mv2gamma <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean) 
}

gamma2mv <- function(shape, rate=NULL, scale=NULL){
  if (is.null(rate)){
    rate <- 1/scale
  }
  
  list(mean = shape/rate,
       var  = shape/rate^2)
}


time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- mv2gamma(mean, var)
    return(rgamma(n, shape = parms$shape, rate = parms$rate))
  } else{
    return(rep(mean, n))
  }
}

time_to_event_lnorm <- function(n, meanlog, sdlog){
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

gen_screening_draws <- function(x){
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}



calc_outcomes <- function(x){
  #browser()
  
  
  # what's the probability of detection at each test time given a value of CT?
  x_ <- x %>%
    mutate(
      ct               = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
      test_p           = case_when(
        assay == "Innova (KCL)"         ~ stats::predict(innova_mod,
                                                         newdata =
                                                           data.frame(ct = ct),
                                                         type =
                                                           "response"
        ),
        assay == "PCR" & ct < 35  ~ 1,
        assay == "PCR" & ct >= 35 ~ 0
      ),
      screen           = runif(n(), 0, 1),
      test_label       = detector(pcr = test_p,  u = screen)
    )
  
  return(x_)
}

calc_outcomes_pre <- function(x){
  #symptomatic self-selection
  #browser()
  x_ <- x %>% 
    mutate(symp_screen = flight_dep_t >= onset_t & flight_dep_t <= end & type=="symptomatic",
           symp_label  = as.logical(symp_screen * rbinom(n=n(),size = 1,prob=synd_screen_prob))) %>% 
    select(-symp_screen)
  
  #pre-board screening?
  x_ <- x_ %>%
    mutate(
      pre_board_t      = ifelse(!is.na(pre_board_assay), flight_dep_t-pre_board_test_delay, NA),
      pre_board_ct     = pmap_dbl(.f = calc_sensitivity, list(model = m, x = pre_board_t)),
      pre_board_p      = case_when(
        pre_board_assay == "Innova (KCL)"         ~ 
          stats::predict(innova_mod,
                         newdata =
                           data.frame(ct = pre_board_ct),
                         type =
                           "response"
          ),
        pre_board_assay == "PCR" &
          pre_board_ct < 35  ~ 1,
        pre_board_assay == "PCR" &
          pre_board_ct >= 35 ~ 0
      ),
      pre_board_screen = runif(n(), 0, 1),
      pre_board_label  = detector(pcr = pre_board_p,  u = pre_board_screen)
    ) 
  
  return(x_)
}


detector <- function(pcr, u = NULL, spec = 1){
  
  if (is.null(u)){
    u <- runif(n = length(pcr))
  }
  
  # true positive if the PCR exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- pcr > u
  
  # false positive if in the top (1-spec) proportion of random draws
  FP <- (pcr == 0)*(runif(n = length(pcr)) > spec)
  
  TP | FP
}


make_delay_label <- function(x,s){
  paste(na.omit(x), s)
}

# for the x axis percentages
trim_percent <- function(x, ...){
  percent(x, ...) %>%
    str_replace(string = ., pattern =  "0+\\%", replacement = "\\%") %>%
    str_replace(string = ., pattern = "\\.\\%", replacement = "\\%")
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_trajectories <- function(n_cases=100, n_sims=100, input, seed=1000,asymp_parms=asymp_fraction){
  
  set.seed(seed)
  #simulate CT trajectories
  
  inf <- data.frame(sim=1:n_sims) %>% 
    mutate(prop_asy = rbeta(n = n(),
                            shape1 = asymp_parms$shape1,
                            shape2 = asymp_parms$shape2)) 
  
  inf %<>%  
    mutate(x = map(.x = prop_asy,
                   .f = ~make_proportions(prop_asy     = .x,
                                          n_cases      = n_cases))) %>% 
    unnest(x) 
  
  
  traj <- inf %>% 
    crossing(start=0) %>% 
    mutate(u = runif(n(),0,1)) %>%
    # duration from: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30172-5/fulltext
    # scaling of asymptomatics taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(end=case_when(type == "symptomatic"  ~ qnormTrunc(p = u, mean=17, 
                                                             sd=approx_sd(15.5,18.6), min = 0),
                         type == "asymptomatic" ~ qnormTrunc(p = u, mean=17*0.6, 
                                                             sd=approx_sd(15.5,18.6), min = 0))) %>% 
    # incubation period from https://bmjopen.bmj.com/content/10/8/e039652.info
    mutate(onset_t=qlnormTrunc(p = u,
                               meanlog=1.63,
                               sdlog=0.5,
                               min = 0,
                               max=end)) %>% 
    pivot_longer(cols = -c(sim,prop_asy,idx,type,u),
                 values_to = "x") %>% 
    # peak CT taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate(y=case_when(name=="start"   ~ 40,
                       name=="end"     ~ 40,
                       name=="onset_t" ~ rnorm(n=n(),mean=22.3,sd=4.2)))
  
  models <- traj %>%
    nest(data = -c(idx,type,u)) %>%  
    dplyr::mutate(
      # Perform loess calculation on each individual 
      m  = purrr::map(data, ~splinefunH(x = .x$x, y = .x$y,
                                        m = c(0,0,0))),
      rx = purrr::map(data, ~range(.x$x)),
      ry = purrr::map(data, ~range(.x$y))) 
  
  models <- models %>% 
    unnest(data) %>%  
    select(-c(y)) %>% 
    pivot_wider(names_from=name,values_from = x) %>%
    mutate(flight_dep_t  = runif(n=n(),min = start,max=end))
  
  return(models)
}

inf_curve_func <- function(m){
  #browser()
  
  x <- data.table(t=seq(0,30,length.out = 51))
  
  #predict CTs for each individual per day
  x$ct <- m(x$t)
  
  #predict culture probability given CTs
  x$culture <-  stats::predict(culture_mod, type = "response", newdata = x)
  
  return(x)
}

inf_threshold <- function(x){
  
  #set 25% infectivity cut-off for infectious period
  x[, `:=`(threshold = 25)][, .(inf_start = min(t[ct < threshold]), 
                                inf_end = max(t[ct < threshold])), 
                            keyby = .(threshold)][, .(inf_start, inf_end)]
}


## just making sure the proportion of cases are secondary or not
make_proportions <- function(prop_asy, n_cases){
  
  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  x <- data.frame(type=rbinom(n=n_cases,size=1,prob = prop_asy)) %>% 
    mutate(type=ifelse(type==1,"asymptomatic","symptomatic"),
           idx=row_number())
  
}


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  gamma2mv(ans$estimate[["shape"]],
           ans$estimate[["rate"]])
}


rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  #browser()
  p_b <- pgamma(q = b, shape = shape, rate = rate)
  p_a <- pgamma(q = a, shape = shape, rate = rate)
  
  u   <- runif(n = n, min = p_a, max = p_b)
  q   <- qgamma(p = u, shape = shape, rate = rate)
  
  return(q)
}


check_unique_values <- function(df, vars){
  # given a data frame and a vector of variables to be used to facet or group, 
  # which ones have length < 1?
  
  l <- lapply(X = vars, 
              FUN =function(x){
                length(unique(df[, x]))
              })
  
  vars[l > 1]
  
}


waning_piecewise_linear <- function(x, ymax, ymin, k, xmax){
  
  if (ymin == ymax){
    Beta = c(0, ymin)
  } else {
    
    Beta <- solve(a = matrix(data = c(xmax, 1,
                                      k,    1),    ncol = 2, byrow = T),
                  b = matrix(data = c(ymin, ymax), ncol = 1))
  }
  
  (x >= 0)*pmin(ymax, pmax(0, Beta[2] + Beta[1]*x))
  
}

waning_points <- function(x, X, Y, log = FALSE){
  
  if (length(X) != length(Y)){
    stop("X and Y must be same length")
  }
  
  if (length(Y) == 1){
    return(rep(Y, length(x)))
  }
  
  if (log){
    Y <- log(Y)
  }
  
  Beta <- solve(a = cbind(X, 1), b = matrix(Y,ncol=1))
  
  Mu <- Beta[2] + Beta[1]*x
  if (log){
    Mu <- exp(Mu)
  }
  (x >= 0)*pmax(0, Mu)
  
}



summarise_simulation <- function(x, faceting, y_labels = NULL){
  
  if(is.null(y_labels)){
    # if none specified, use all.
    y_labels_names <- grep(x=names(x), pattern="^trans_pot_", value = T)
  } else {
    y_labels_names <- names(y_labels)
  }
  
  all_grouping_vars <- all.vars(faceting)
  
  # if (!any(grepl(pattern = "type", x = all_grouping_vars))){
  #   all_grouping_vars <- c(all_grouping_vars, "type")
  # }
  
  x_summaries <-
    as.list(y_labels_names) %>%
    set_names(., .) %>%
    lapply(X = ., 
           FUN = function(y){
             make_quantiles(x,
                            y_var = y, 
                            vars = all_grouping_vars)})
  
  if (any(grepl(pattern = "type", x = all_grouping_vars))){
    
    x_summaries_all <- as.list(y_labels_names) %>%
      set_names(., .) %>%
      lapply(X = ., 
             FUN = function(y){
               make_quantiles(
                 mutate(x,
                        type = "all"),
                 y_var = y, 
                 vars = all_grouping_vars)})
    
    x_summaries <- map2(.x = x_summaries,
                        .y = x_summaries_all,
                        .f = ~bind_rows(.x, .y))
    
  }
  
  bind_rows(x_summaries, .id = "yvar")
  
}

calc_sensitivity <- function(model, x){
  #browser()
  if(!is.na(x)){
    s <- model(x)
  } else {
    s <- NA_real_
  }
  
  return(s)
}



read_results <- function(results_path){
  #browser()
  list(here::here("results", results_path, "results.rds"),
       here::here("results", results_path, "input.rds")) %>%
    map(read_rds) %>%
    map(bind_rows) %>%
    {inner_join(.[[1]], .[[2]])}
}


test_times <- function(multiple_tests,tests,flight_dep_t,flight_arr_t,quar_dur,sampling_freq = 1, n_tests){
  #browser()
  
  if(multiple_tests){
    
    test_timings <- data.table(test_t = seq(from=flight_arr_t,to=(flight_arr_t+n_tests-1),by=sampling_freq)) 
    
    test_timings[, `:=`(test_no = paste0("test_", seq_len(.N)))]
    
  } else {
    test_timings <- data.table(test_t=flight_arr_t+quar_dur)
    
    test_timings[, `:=`(test_no = paste0("test_",  seq_len(.N)))]
  }
  
  return(test_timings)
}

earliest_pos <- function(df){
  #
  #browser()
  df <- as.data.table(df)
  
  df <- df[test_label==TRUE]
  
  if (nrow(df) == 0L){
    data.table(test_no="None",test_p=0,test_t=Inf)
  } else {
    df[, .(test_no, test_p, test_t)][, .SD[order(test_t)][frankv(test_t, ties.method = "min", na.last = "keep") <= 1L]]
  }
}

calc_inf_period <- function(x){
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), x %>% slice(1) %>%  pull(scenario)))
  x_ <- as.data.table(x)
  
  x_[, inf_curve := map(.f=inf_curve_func,.x=m)]
  x_[, inf_period := map(.f=inf_threshold,.x=inf_curve)]
  
  x_ <- unnest(x_,inf_period)
  
}


calc_auc <- function(x){
  #browser()
  
  x_ <- as.data.table(x)
  
  x_[, inf_curve := map(.f=inf_curve_func,.x=m)]
  x_[, inf_period := map(.f=inf_threshold,.x=inf_curve)]
  
  x_[, `:=` (adhering_quar = rbinom(n = .N, size = 1,  prob = adherence_quar), 
             adhering_symp = rbinom(n = .N, size = 1, prob = adherence_symp), 
             adhering_test = rbinom(n = .N, size = 1, prob = adherence_test))]
  
  x_ <- unnest(x_,inf_period)
  
  x_ <- x_ %>% 
    rowwise() %>% 
    mutate(
      pre_symp_auc       = case_when(symp_label~
                                       MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                                                 from = flight_dep_t),
                                     TRUE ~ 0),
      pre_board_test_auc = case_when(pre_board_label~
                                       MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                                                 from = flight_dep_t),
                                     TRUE ~ 0),
      symp_auc           = case_when(symp_label ~ 0,
                                     pre_board_label~0,
                                     type == "asymptomatic" ~0,
                                     type== "symptomatic" ~
                                       MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                                                 from = pmax(flight_dep_t,onset_t),
                                                 to   = pmax(flight_dep_t,symp_end_t)) * adhering_symp
      ),
      quar_auc           = case_when(
        quar_dur==0~0,
        multiple_tests ~ 0,
        symp_label ~ 0,
        pre_board_label~0,
        #If symptomatic before quarantine, no quarantine
        !multiple_tests & onset_t < flight_arr_t & type=="symptomatic" ~0,
        #If symptomatic during quarantine, quarantine ends at onset
        !multiple_tests & onset_t> flight_arr_t & onset_t < quar_end_t & type=="symptomatic" ~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = flight_arr_t,
                    to   = onset_t) * adhering_quar,
        #If symptomatic after quarantine, quarantine ends at standard time
        !multiple_tests & onset_t>quar_end_t&type=="symptomatic"~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = flight_arr_t,
                    to   = quar_end_t) * adhering_quar,
        #If asymptomatic, quarantine is completed in full
        !multiple_tests&type=="asymptomatic" ~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = flight_arr_t,
                    to   = quar_end_t) * adhering_quar),
      test_auc           = case_when(
        !tests ~ 0,
        symp_label ~ 0,
        pre_board_label~0,
        #If testing and symptoms occur before the test, no-post test isolation
        tests & onset_t < earliest_t &type =="symptomatic" ~ 0,
        #If testing and symptoms occur after a positive test, end post-positive test self-isolation
        tests & onset_t > earliest_t & onset_t < test_iso_end_t & type == "symptomatic" ~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = earliest_t,
                    to   = onset_t) * adhering_test,
        #If testing and symptoms occur after the end of post-positive test isolation, complete full self-isolation
        tests & onset_t > earliest_t & onset_t > test_iso_end_t & type == "symptomatic"~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = earliest_t,
                    to   = test_iso_end_t) * adhering_test,
        #If testing and asymptomatic, complete full self-isolation
        tests & type == "asymptomatic"~
          MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                    from = earliest_t,
                    to   = test_iso_end_t) * adhering_test),
      max_auc = pre_board_test_auc + pre_symp_auc + symp_auc + quar_auc + test_auc,
      tot_auc = MESS::auc(x=inf_curve$t,y=inf_curve$culture,
                          from = pmax(flight_dep_t,0))
    ) 
  
}

approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
}

arrivals_fun <- function(x, trav_vol, sims=1){
  
  weekly_arrivals <- rbinom(n = sims, size = trav_vol, prob = 7/30)  
  
  infected_arrivals <- rmultinom(n = sims, 
                                 size = weekly_arrivals,
                                 prob = x$prev) %>%
    unlist %>%
    {data.frame(arrivals = ., type = x$type)} %>%
    filter(type != "uninfected")
  
  #rbinom(n = sims, prob=prev, size = weekly_arrivals)
  
  infected_arrivals
}

make_col_pal <- function(my_palette = grey.colors, x, ...){
  
  #browser()

  
  # ensure that all characters have levels as factors
  if ("character" %in% class(x)){
    x <- factor(x)
  }
  
  # turn numbers into factors, respecting order not alphabetical
  if ("numeric" %in% class(x)){
    x <- factor(as.character(x), levels = sort(unique(x)), ordered = T)
  }
  
  # how many colours and how are they labelled
  n         <- nlevels(x)
  col_names <-  levels(x)
  
  if ("character" %in% class(my_palette)){
    # do we already have a valid colour palette?
    if(all(my_palette %in% colors() | 
           grepl(pattern = "#[0-9A-F]{6}", x = my_palette))){
      col_pal <- my_palette[col_names]
    }
    
    # have we been passed a named RColorBrewer palette?
    if(my_palette[[1]] %in% rownames(RColorBrewer::brewer.pal.info)){
      if (length(my_palette) > 1L){
        warning("Using first entry in colour palette")
      }
      col_pal <- RColorBrewer::brewer.pal(n    = n,
                                          name = my_palette[[1]])
    }
    
  }
  
  # have we been passed a colour palette generating function?
  if ("function" %in% class(my_palette)){
    col_pal <- my_palette(n, ...)
  }
  
  col_pal
  
}




rr_func <- function(x=results_df,
                    baseline,
                    x_var=NULL,
                    row_vars=NULL,
                    col_vars=NULL,
                    group_var=stringency,
                    probs = c(0.025,0.25,0.5,0.75,0.975),
                    log=T, my_palette = 'Dark2'){
  #browser()
  
  sim_dots <- sym("sim")
  x_dots <-  enquo(x_var) 
  row_dots <- enquo(row_vars)
  col_dots  <- enquo(col_vars)
  group_dots <- enquo(group_var)
  dots  <- enquos(sim_dots,x_var,row_vars,col_vars,group_var)
  
  
  x_ <- x %>%
    filter(!is.na({{x_var}})) %>%
    group_by(!!!dots) %>% 
    #filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
    summarise(prop=1-(sum(max_auc)/sum(tot_auc))) %>% 
    left_join(baseline,by="sim") %>% 
    mutate(prop_ratio=prop/baseline_prop) %>% 
    replace_na(list(prop_ratio=1)) 
  
  col_pal <- RColorBrewer::brewer.pal(n = pull(x_, !!group_dots) %>% 
                                        nlevels(),
                                      name = my_palette)
  
  names(col_pal) <- x_ %>% 
    pull({{group_var}}) %>% 
    levels()
  
  
  x_plot <- x_ %>%
    group_by(!{{x_var}},!!!row_dots,!!!col_dots,!{{group_var}}) %>% 
    nest() %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .$prop_ratio,
                                                probs = probs)),
           Mean = map_dbl(.x=data,
                          ~mean(.$prop_ratio,na.rm=T)),
           SD   = map_dbl(.x=data,
                          ~sd(.$prop_ratio,na.rm = T))) %>%
    unnest_wider(Q) %>% 
    dplyr::select(-data) %>% 
    ggplot(aes(x = factor({{x_var}}), y = `50%`)) + 
    geom_hline(aes(yintercept=1),colour="grey",linetype="dashed")+
    geom_linerange(aes(ymin = `2.5%`,
                       ymax = `97.5%`,
                       colour={{group_var}}),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
    geom_linerange(aes(ymin = `25%`,
                       ymax = `75%`,
                       colour={{group_var}}),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
    geom_point(aes(y = `50%`,colour={{group_var}}),
               #pch="-",
               size=1.5,
               position=position_dodge(width=0.5)) +
    facet_nested(nest_line=T,
                 rows = vars(!!row_dots), 
                 cols = vars(!!col_dots),
                 labeller = labeller(
                   type = capitalize,
                   sens_LFA=sens_scaling_labels,
                   delay_scaling = delay_scaling_labeller,
                   pre_board_test = pre_board_labels,
                   adherence =
                     c("full" = "Full adherence",
                       "literature" = "Adherence from literature")
                 )) +
    scale_colour_manual(breaks = names(col_pal),
                        values= col_pal,
                        name=NULL)+
    plotting_theme
  
  x_plot <- x_plot + labs(x=case_when(quo_text(x_dots)=="quar_dur"~"Quarantine duration (days)",
                                      TRUE ~ "Number of days of testing"),
                          y=case_when(baseline$name=="baseline_high"~"Transmission averted compared tosyndromic\nscreening and 14 day quarantine on arrival",
                                      TRUE~"Transmission averted compared to syndromic\nscreening and symptomatic self-isolation only"))
  
  x_plot <- x_plot + scale_y_continuous(trans = case_when(log==TRUE~"log10",
                                                          TRUE~"identity"),
                                        breaks = breaks_width(0.1))
  
  x_summaries <- x_ %>%  
    group_by(!{{x_var}},!!!row_dots,!!!col_dots,!{{group_var}}) %>% 
    nest() %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .$prop_ratio,
                                                probs = probs)),
           Mean = map_dbl(.x=data,
                          ~mean(.$prop_ratio,na.rm=T)),
           SD   = map_dbl(.x=data,
                          ~sd(.$prop_ratio,na.rm = T))) %>%
    unnest_wider(Q) %>% 
    mutate_at(.vars = vars(contains("%")), .funs = txtRound,digits=2) %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
    select(`50%`,ui)
  
  return(list(summaries=x_summaries,plot=x_plot))
  
}


plotting_func <- function(x=results_df,
                          x_var=NULL,
                          group_var=assay,
                          faceting=NULL,
                          probs = c(0.025,0.25,0.5,0.75,0.975),
                          include_symp=F,
                          my_palette = 'Dark2',
                          zero_line = T){
 
  if(!is.null(faceting)){
    col_facets <- formula.tools::rhs.vars(faceting)
    row_facets <- formula.tools::lhs.vars(faceting)
  } else {
    col_facets <- NULL
    row_facets <- NULL
  }
  
  if(include_symp){
    x_ <- x %>%
      filter(!is.na({{ x_var }})) %>%
      group_by_at(.vars=vars(sim,{{x_var}},
                             {{group_var}},
                             all_of(col_facets),
                             all_of(row_facets))) %>% 
      summarise(prop=-sum(max_auc,na.rm = T)/sum(tot_auc,na.rm = T))
  } else {
    x_ <- x %>%
      filter(!is.na({{ x_var }})) %>%
      group_by_at(.vars=vars(sim,{{x_var}},
                             {{group_var}},
                             all_of(col_facets),
                             all_of(row_facets))) %>% 
      summarise(prop=-sum(max_auc-symp_auc,na.rm = T)/sum(tot_auc-symp_auc,na.rm = T))
  }
  
  #browser()
  
  pal_vals <- pull(distinct(ungroup(x_), {{group_var}}), {{group_var}})
  
  col_pal <- make_col_pal(my_palette = my_palette, x = pal_vals)
  
  x_plot <- x_ %>%
    group_by_at(.vars=vars({{x_var}},{{group_var}},all_of(col_facets),all_of(row_facets))) %>% 
    nest() %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .$prop,
                                                probs = probs)),
           Mean = map_dbl(.x=data,
                          ~mean(.$prop,na.rm=T)),
           SD   = map_dbl(.x=data,
                          ~sd(.$prop,na.rm = T))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    ggplot(aes(x = factor({{x_var}}), y = `50%`)) + 
    geom_linerange(aes(ymin = `2.5%`,
                       ymax = `97.5%`,
                       colour={{group_var}}),
                   position=position_dodge(width=0.5),
                   size=1,
                   alpha=0.5) +
    geom_linerange(aes(ymin = `25%`,
                       ymax = `75%`,
                       colour={{group_var}}),position=position_dodge(width=0.5),size=1.2,alpha=0.7)+
    geom_point(aes(y = `50%`,colour={{group_var}}),
               shape=21,
               fill="white",
               size=1,
               position=position_dodge(width=0.5)) +
    plotting_theme+
    scale_colour_manual(values = col_pal,
                        breaks = names(col_pal),
                        name=NULL)
  
  if(!is.null(faceting)){
    x_plot <- x_plot+
    facet_nested(nest_line=T,
                 faceting,
                 labeller = labeller(
                   type           = capitalize,
                   sens_LFA       = sens_scaling_labels,
                   delay_scaling  = delay_scaling_labeller,
                   pre_board_test = pre_board_labels,
                   pre_board_assay = function(x){str_wrap(string = x,
                                                          width = 15)},
                   adherence      =
                     c("full" = "Full adherence",
                       "literature" = "Adherence values\nfrom literature")
                 )) 
  }
  
  if (zero_line) {
    x_plot <- x_plot + geom_hline(aes(yintercept=0),linetype="dashed",colour="grey")
  }
  
  xlab_name <- names(select(ungroup(x_), {{x_var}}))
  
  x_plot <- x_plot + labs(
    x=case_when(
      xlab_name == "quar_dur" ~ "Quarantine duration (days)",
      xlab_name == "n_tests"  ~ "Number of days of testing",
      TRUE ~ xlab_name),
    y="Transmission averted")
  
  x_summaries <- x_ %>%  
    group_by_at(.vars=vars({{x_var}},{{group_var}},all_of(col_facets),all_of(row_facets))) %>% 
    nest() %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .$prop,
                                                probs = probs)),
           Mean = map_dbl(.x=data,
                          ~mean(.$prop,na.rm=T)),
           SD   = map_dbl(.x=data,
                          ~sd(.$prop,na.rm = T))) %>%
    unnest_wider(Q) %>% 
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(ui, c(`97.5%`,`2.5%`), sep= ", ") %>% 
    mutate(ui=paste0("(95% UI: ",ui,")")) %>% 
    select(`50%`,ui)
  
  return(list(summaries=x_summaries,plot=x_plot))
  
}

gamma_draw <- function(low,high){
  #browser()
  if(is.na(low)|is.na(high)){
    draw <- NA_real_
  }else{
    x_ <- gamma.parms.from.quantiles(q=c(low,high))
    
    draw <- rgamma(n=1,shape=x_$shape,rate=x_$rate)
  }
  return(draw)
  
}

`%!in%` <- negate(`%in%`)


gamma.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                       precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1.1, 1.1), plot=F, plot.xlim=numeric(0))
{
  # Version 1.0.2 (December 2012)
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Bélisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/GammaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dgamma(x, shape=theta[1], rate=theta[2])}
  F.inv <- function(x, theta){qgamma(x, shape=theta[1], rate=theta[2])}
  f.cum <- function(x, theta){pgamma(x, shape=theta[1], rate=theta[2])}
  f.mode <- function(theta){shape <- theta[1]; rate <- theta[2]; mode <- ifelse(shape>1,(shape-1)/rate,NA); mode}
  theta.from.moments <- function(m, v){shape <- m*m/v; rate <- m/v; c(shape, rate)}
  dens.label <- 'dgamma'
  parms.names <- c('shape', 'rate')
  
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        k <- min(last.theta/change)
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(shape=parms$theta[1], rate=parms$theta[2], scale=1/parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}

beta.parms.from.quantiles <- function(q, p=c(0.025,0.975),
                                      precision=0.001, derivative.epsilon=1e-3, start.with.normal.approx=T, start=c(1, 1), plot=F)
{
  # Version 1.3 (February 2017)
  #
  # Function developed by 
  # Lawrence Joseph and pbelisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
  #
  # Please refer to our webpage for details on each argument.
  
  f <- function(x, theta){dbeta(x, shape1=theta[1], shape2=theta[2])}
  F.inv <- function(x, theta){qbeta(x, shape1=theta[1], shape2=theta[2])}
  f.cum <- function(x, theta){pbeta(x, shape1=theta[1], shape2=theta[2])}
  f.mode <- function(theta){a <- theta[1]; b <- theta[2]; mode <- ifelse(a>1, (a-1)/(a+b-2), NA); mode}
  theta.from.moments <- function(m, v){a <- m*m*(1-m)/v-m; b <- a*(1/m-1); c(a, b)}
  plot.xlim <- c(0, 1)
  
  dens.label <- 'dbeta'
  parms.names <- c('a', 'b')
  
  if (length(p) != 2) stop("Vector of probabilities p must be of length 2.")
  if (length(q) != 2) stop("Vector of quantiles q must be of length 2.")
  p <- sort(p); q <- sort(q)
  
  #_____________________________________________________________________________________________________
  
  print.area.text <- function(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim, M=30, M0=50)
  {
    par.usr <- par('usr')
    par.din <- par('din')
    
    p.string <- as.character(round(c(0,1) + c(1,-1)*p.check, digits=4))
    str.width <- strwidth(p.string, cex=cex)
    str.height <- strheight("0", cex=cex)
    
    J <- matrix(1, nrow=M0, ncol=1)
    
    x.units.1in <- diff(par.usr[c(1,2)])/par.din[1]
    y.units.1in <- diff(par.usr[c(3,4)])/par.din[2]
    aspect.ratio <- y.units.1in/x.units.1in
    
    # --- left area  -----------------------------------------------------------
    
    scatter.xlim <- c(max(plot.xlim[1], par.usr[1]), q[1])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=0, to=p[1], length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-1]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the right of the mode, if any
    w <- which(x>mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[1]+str.width[1]) <= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y, fromLast=T))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[1], adj=c(1,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[1], mean(par.usr[c(3,4)]), labels=p.string[1], col='gray', cex=cex, srt=90, adj=c(1,0))
    }
    
    # --- right area  ----------------------------------------------------------
    
    scatter.xlim <- c(q[2], plot.xlim[2])
    scatter.ylim <- c(0, par.usr[4])
    x <- seq(from=scatter.xlim[1], to=scatter.xlim[2], length=M)
    y <- seq(from=scatter.ylim[1], to=scatter.ylim[2], length=M)
    x.grid.index <- rep(seq(M), M)
    y.grid.index <- rep(seq(M), rep(M, M))
    grid.df <- f(x, theta)
    
    # Estimate mass center
    tmp.p <- seq(from=p[2], to=f.cum(plot.xlim[2], theta), length=M0)
    tmp.x <- F.inv(tmp.p, theta)
    h <- f(tmp.x, theta)
    mass.center <- c(mean(tmp.x), sum(h[-length(h)]*diff(tmp.x))/diff(range(tmp.x)))
    
    # Identify points under the curve
    # (to eliminate them from the list of candidates)
    gridpoint.under.the.curve <- y[y.grid.index] <= grid.df[x.grid.index]
    w <- which(gridpoint.under.the.curve)
    x <- x[x.grid.index]; y <- y[y.grid.index]
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points to the left of the mode, if any
    w <- which(x<mode)
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # Eliminate points for which the text would fall out of the plot area
    w <- which((par.usr[2]-str.width[2]) >= x & (y + str.height) <= par.usr[4])
    x <- x[w]; y <- y[w]
    
    # For each height, eliminate the closest point to the curve
    # (we want to stay away from the curve to preserve readability)
    w <- which(!duplicated(y))
    if (length(w)){x <- x[-w]; y <- y[-w]}
    
    # For each point, compute distance from mass center and pick the closest point
    d <- ((x-mass.center[1])^2) + ((y-mass.center[2])/aspect.ratio)^2
    w <- which.min(d)
    x <- x[w]; y <- y[w]
    
    if (length(x))
    {
      text(x, y, labels=p.string[2], adj=c(0,0), col='gray', cex=cex)
    }
    else
    {
      text(plot.xlim[2], mean(par.usr[c(3,4)]), labels=p.string[2], col='gray', cex=cex, srt=-90, adj=c(1,0))
    }
  }
  
  # ......................................................................................................................................
  
  Newton.Raphson <- function(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start)
  {
    Hessian <- matrix(NA, 2, 2)
    
    if (start.with.normal.approx)
    {
      # Probably not a very good universal choice, but proved good in most cases in practice
      m <-  diff(q)/diff(p)*(0.5-p[1]) + q[1]
      v <- (diff(q)/diff(qnorm(p)))^2
      theta <- theta.from.moments(m, v)
    }
    else theta <- start
    
    
    change <- precision + 1
    niter <- 0
    # Newton-Raphson multivariate algorithm
    while (max(abs(change)) > precision)
    {
      Hessian[,1] <- (f.cum(q, theta) - f.cum(q, theta - c(derivative.epsilon, 0))) / derivative.epsilon
      Hessian[,2] <- (f.cum(q, theta) - f.cum(q, theta - c(0, derivative.epsilon))) / derivative.epsilon
      
      f <- f.cum(q, theta) - p
      change <- solve(Hessian) %*% f
      last.theta <- theta
      theta <- last.theta - change
      
      # If we step out of limits, reduce change
      
      if (any(theta<0))
      {
        w <- which(theta<0)
        k <- min(last.theta[w]/change[w])
        theta <- last.theta - k/2*change
      }
      
      niter <- niter + 1
    }
    
    list(theta=as.vector(theta), niter=niter, last.change=as.vector(change))
  }
  
  # ...............................................................................................................
  
  plot.density <- function(p, q, f, f.cum, F.inv, mode, theta, plot.xlim, dens.label, parms.names, cex)
  {
    if (length(plot.xlim) == 0)
    {
      plot.xlim <- F.inv(c(0, 1), theta)
      
      if (is.infinite(plot.xlim[1]))
      {
        tmp <- min(c(0.001, p[1]/10))
        plot.xlim[1] <- F.inv(tmp, theta)
      }  
      
      if (is.infinite(plot.xlim[2]))
      {
        tmp <- max(c(0.999, 1 - (1-p[2])/10))
        plot.xlim[2] <- F.inv(tmp, theta)
      }
    }
    plot.xlim <- sort(plot.xlim)
    
    
    x <- seq(from=min(plot.xlim), to=max(plot.xlim), length=1000)
    h <- f(x, theta)
    x0 <- x; f0 <- h
    ylab <- paste(c(dens.label, '(x, ', parms.names[1], ' = ', round(theta[1], digits=5), ', ', parms.names[2], ' = ', round(theta[2], digits=5), ')'), collapse='')
    plot(x, h, type='l', ylab=ylab)
    
    # fill in area on the left side of the distribution
    x <- seq(from=plot.xlim[1], to=q[1], length=1000)
    y <- f(x, theta)
    x <- c(x, q[1], plot.xlim[1]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # fill in area on the right side of the distribution
    x <- seq(from=max(plot.xlim), to=q[2], length=1000)
    y <- f(x, theta)
    x <- c(x, q[2], plot.xlim[2]); y <- c(y, 0, 0)
    polygon(x, y, col='lightgrey', border='lightgray')
    # draw distrn again
    points(x0, f0, type='l')
    h <- f(q, theta)
    points(rep(q[1], 2), c(0, h[1]), type='l', col='orange')
    points(rep(q[2], 2), c(0, h[2]), type='l', col='orange')
    # place text on both ends areas
    print.area.text(p, p.check, q, f, f.cum, F.inv, theta, mode, cex, plot.xlim)  
    
    xaxp <- par("xaxp")
    x.ticks <- seq(from=xaxp[1], to=xaxp[2], length=xaxp[3]+1)
    q2print <- as.double(setdiff(as.character(q), as.character(x.ticks)))
    
    mtext(q2print, side=1, col='orange', at=q2print, cex=0.6, line=2.1)
    points(q, rep(par('usr')[3]+0.15*par('cxy')[2], 2), pch=17, col='orange')
  }
  
  #________________________________________________________________________________________________________________
  
  
  parms <- Newton.Raphson(derivative.epsilon, precision, f.cum, p, q, theta.from.moments, start.with.normal.approx, start=start)
  p.check <- f.cum(q, parms$theta)
  
  if (plot) plot.density(p, q, f, f.cum, F.inv, f.mode(parms$theta), parms$theta, plot.xlim, dens.label, parms.names, 0.8)
  
  list(a=parms$theta[1], b=parms$theta[2], last.change=parms$last.change, niter=parms$niter, q=q, p=p, p.check=p.check)
}


library(rlang)
unnest_dt <- function(tbl, col) {
  tbl <- as.data.table(tbl)
  col <- ensyms(col)
  clnms <- syms(setdiff(colnames(tbl), as.character(col)))
  tbl <- as.data.table(tbl)
  tbl <- eval(
    expr(tbl[, as.character(unlist(!!!col)), by = list(!!!clnms)])
  )
  colnames(tbl) <- c(as.character(clnms), as.character(col))
  tbl
}


label_template <- function(x, s = ""){
  paste(s,
        percent(parse_number(x),
                accuracy = 1))
}

`%!in%` <- Negate(`%in%`)
