## functions used for ploting

# colours
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
covid_pal2 <- set_names(covid_pal, c("All", "Asymptomatic", "Symptomatic"))
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

# assay_colours <- RColorBrewer::brewer.pal(n = length(assay_labeller),
#                                          name = "Dark2") %>%
#   set_names(names(assay_labeller))

extrafont::loadfonts()
pdf.options(useDingbats=FALSE)

plotting_theme <- 
  theme_minimal(base_family = "Lato") +
  theme(panel.border    = element_rect(fill=NA),
        legend.position = "bottom",
        legend.box      = "horiztonal",
        legend.text.align = 0,
        axis.ticks = element_line(),
        panel.grid.minor = element_blank())

infectivity_labels <-
  c("infectivity_post" =
      "Transmission potential of secondary cases \nafter release",
    "infectivity_averted" = 
      "Transmission potential of secondary cases \naverted as a result of quarantine and testing",
    "infectivity_avertable" = 
      "Remaining transmission potential of secondary cases \n post-tracing",
    "infectivity_quar" = 
      "Transmission potential in community\ndue to imperfect quarantine adherence",
    "infectivity_pre" =
      "Transmission potential of secondary cases \nprior to being traced",
    "infectivity_total" = 
      "Transmission potential of secondary cases \nin community compared to no quarantine or testing",
    "infectivity_mass"  = 
      "Truncated transmission potential\ndistribution mass",
    "infectivity_avertable" = 
      "Transmission potential occurring\nafter quarantine starts",
    "infectivity_post_release_onset" =
      "Transmission potential averted during\npost-quarantine self-isolation"
  )



test_labeller <- function(x){
  
  mutate(x,
         stringency = factor(stringency,
                             levels = c("none",
                                        "one",
                                        "two"),
                             labels = c("No test",
                                        "One test",
                                        "Two tests"),
                             ordered = T))
}

type_labels <- c("asymptomatic" =
                   "Asymptomatic",
                 "symptomatic" =
                   "Pre-symptomatic")

index_test_labeller <- function(x, newline = FALSE){
  paste0("Delay between index case's\nonset and having a test:",
         ifelse(newline, "\n", " "),
         x,
         ifelse(x == 1, " day", " days"))
}

# delay_scaling_labeller <- function(x, newline = FALSE){
#   paste0("Contact tracing delay",
#          ifelse(newline, "\n", " "),
#          "scaling factor: ",
#          x)
# }
delay_scaling_labeller <- function(x){
  dplyr::case_when(
    x==0   ~ "Instant T&T (0 days)",  
    x==0.5 ~ "T&T delays halved (1.5 days)",
    x==1   ~ "Observed T&T delays (3 days)",
    TRUE   ~ "Unknown")
}

sens_scaling_labels <- c(
  "higher" = "Oxford/PHE evaluation",
  "lower"  = "Liverpool pilot (baseline)")

pre_board_labels <- function(x){
  case_when(x==TRUE~"Pre-board testing",
            x==FALSE~"No pre-board testing")
}


waning_labeller <- function(x){
  #paste("Adherence to quarantine guidance:\n",
  dplyr::case_when( 
    x == "waning_none"             ~ "Constant 100% adherence",
    x == "waning_constant"         ~ "Constant 75% adherence",
    x == "waning_canada_total"     ~ "Exponential waning from 100%",
    x == "waning_canada_community" ~ "Exponential decay (community only)",
    TRUE ~ "Unknown")
  #)
}

percentage <- function(x, ...){
  
  ans <- scales::percent(x, ...)
  ans[x > 1] <- ""
  
  ans
  
}

pretty_percentage <- function(x){
  ans <- pretty(x)
  ans[ans <= 1]
}

make_release_figure <- function(x_summaries,
                                input,
                                xlab = "Days in quarantine",
                                ylab = "",
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                hline = 0,
                                faceting = NULL,
                                percent = FALSE){
  
  x_summaries %<>% test_labeller # should this be in the facet call?
  
  # how to do presymptomatic
  
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x_summaries %<>% mutate(type = factor(type,
                                          levels = c("asymptomatic",
                                                     "symptomatic"),
                                          labels = c("Asymptomatic",
                                                     "Presymptomatic")))
  }
  
  
  figure <-  
    ggplot(data=x_summaries, aes(x = time_since_exp, 
                                 y = `50%`, 
                                 color = stringency)) +
    geom_hline(aes(yintercept=1), linetype=hline)+
    geom_linerange(aes(ymin  = `2.5%`, 
                       ymax  = `97.5%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.3,
                   size      = 3) +
    geom_linerange(aes(ymin  = `25%`,
                       ymax  = `75%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.5,
                   size      = 3) +
    geom_point(pch           = "-",
               size          = 12,
               position      = position_dodge2(width = 1),
               aes(y         = `50%`,
                   group     = stringency)) +
    scale_x_continuous(breaks = breaks_width(2))+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme_minimal()+
    theme(axis.ticks = element_line(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.position = "bottom",
          strip.placement = "outside") +
    ggplot2::labs(x = xlab,
                  y = ylab) +
    xlab("Days in quarantine\n(including 1 day delay on testing results)")
  
  figure <- figure + 
    facet_nested(nest_line = TRUE,
                 facets = faceting,
                 labeller = labeller(index_test_delay = index_test_labeller,
                                     delay_scaling    = delay_scaling_labeller,
                                     waning           = waning_labeller))
  
  
  return(figure)
  
}

plot_data <- function(input, 
                      x_summaries,
                      main_scenarios = NULL){
  
  dat <- x_summaries  %>%
    inner_join(input) %>% # should really carry these through when summarising
    mutate(time_since_exp_ = 
             ifelse(is.na(time_since_exp),
                    0,
                    time_since_exp),
           time_in_iso = 
             first_test_delay + 
             time_since_exp_+
             screening)
  
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% dplyr::select(-one_of("released_test")) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    #tidyr::nest(data = -c(first_test_delay, time_since_exp)) %>%
    dplyr::mutate(delays = paste(first_test_delay, "&",
                                 first_test_delay + time_since_exp)) %>%
    #tidyr::unnest(data) %>%
    dplyr::mutate(time_in_iso = factor(time_in_iso, 
                                       levels = sort(unique(.$time_in_iso)),
                                       ordered = T)) %>%
    dplyr::filter(M!=0) %>%  # if the mean is zero, this group is empty
    return
}

make_arrivals_table <- function(x, table_vars = c("country")){
  x %>%
    mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
    group_by_at(.vars = vars(sim, one_of(table_vars))) %>%
    count(.drop = F) %>%
    group_by_at(.vars = vars(-n, -sim)) %>%
    nest() %>%
    mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    group_by_at(.vars = vars(one_of(table_vars))) %>%
    transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
    mutate_at(.vars = vars(-c(country, value)), 
              .funs = str_to_title) %>%
    spread(country, value)
}


save_plot <- function(plot   = ggplot2::last_plot(),
                      prefix = stringi::stri_rand_strings(1, length = 8),
                      base   = NULL, # to help identify which figure in the paper
                      device = NULL, # additional devices, e.g. "png" 
                      width  = 210, 
                      height = 210,
                      dpi    = 300,
                      units  = "mm"){
  #browser()
  file <- paste0("results/",
                 paste(prefix, base, sep = "_"),
                 ".pdf")
  
  #if(Sys.info()["sysname"]=="Windows"){
  cairo_pdf(filename = file,
        width  = measurements::conv_unit(x = width,  units, "inch"),
        height = measurements::conv_unit(x = height, units, "inch"))
  # }else{
  # pdf(file = file,
  #     width  = measurements::conv_unit(x = width,  units, "inch"),
  #     height = measurements::conv_unit(x = height, units, "inch"),
  #     useDingbats = FALSE
  #     )
  # }
  
  print(plot)
  
  dev.off()
  
  if (length(device) > 0L){
    #img <- pdftools::pdf_render_page(file, dpi = dpi)
    
    img <- magick::image_read_pdf(file, density = dpi)
    
    purrr::map(.x = device,
               .f = 
                 ~magick::image_convert(image = img, format = .x) %>%
                 magick::image_write(., path = sub(pattern = "pdf",
                                                   replacement = .x,
                                                   x = file)))
  }
}



make_days_plots <- 
  function(x, 
           main_scenarios   = NULL,
           plot             = TRUE,
           log_scale        = FALSE,
           text_size        = 2.5,
           xlab             = "Days since exposure\n(including 1 day delay on testing results)",
           sum              = FALSE,
           y_labels         = NULL, # pass in y_vars as a named list
           faceting         = NULL,
           dir              = stringi::stri_rand_strings(1, 8),
           base             = stringi::stri_rand_strings(1, 8)){
    
    if (!dir.exists(paste0("results/",dir))){
      dir.create(paste0("results/",dir))
    }
    
    all_grouping_vars <- all.vars(faceting)
    
    if (sum){
      y_labels <- sub("^Average", "Total", y_labels)
    } 
    
    x_days_summaries <-
      as.list(names(y_labels)) %>%
      lapply(X = ., FUN = function(y){
        map_df(.x = x,
               ~make_released_time_quantiles(.x,
                                             y_var = y, sum = sum,
                                             vars = all_grouping_vars))})
    ## end summaries
    
    fig_data <- x_days_summaries %>% 
      map(~plot_data(input = input, # should we still pass it in? recover?
                     x_summaries = 
                       .x,
                     main_scenarios))
    if (plot){
      
      figs <- map2(
        .x = fig_data,
        .y = y_labels,
        .f = ~make_release_figure(
          x         = .x,
          #input     = input,
          xlab      = xlab,
          text_size = text_size,
          ylab      = .y,
          faceting  = faceting,
          percent   = TRUE) )
      
      
      if (length(figs) > 1){
        legend <- cowplot::get_legend(figs[[1]])
        figs   <- lapply(X = figs, function(x){x + theme(legend.position = "none")})
        
        fig    <- cowplot::plot_grid(plotlist = figs, nrow = 1,
                                     labels = LETTERS[1:length(figs)])
        
        fig    <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .2))
        
      } else {
        fig <- figs[[1]]
      }
      
      
      list("png", "pdf") %>%
        map(~ggsave(filename = paste0("results/",dir,"/days_plots_",base,".",.x),
                    plot=fig,
                    width  = 60*nrow(distinct(fig_data[[1]][,get.vars(rhs(faceting))]))*
                      length(fig_data), 
                    height = 120*nrow(distinct(fig_data[[1]][,get.vars(lhs(faceting))])), 
                    dpi = 600,
                    units = "mm",
                    device = ifelse(.x == "pdf",
                                    cairo_pdf,
                                    "png")))
    }
    
    return(
      set_names(fig_data, sub(pattern = "infectivity_", 
                              replacement = "", x = names(y_labels)))
      
    )
    
  }


summarise_results <- function(x, reduction = TRUE){
  if (!is.logical(reduction)){
    stop("reduction must be logical")
  }
  x <- mutate_at(x, 
                 .vars = vars(contains("%")),
                 .funs = function(x){
                   reduction*(1 - x) +
                     (1 - reduction)*x})
  if (reduction){
    percentages <- grep(x = names(x), pattern = "%")
    x[,rev(percentages)] <- x[,percentages]  
  }
  
  x
  
}

show_results <- function(x, reduction = TRUE){
  dplyr::select(x, delays,
                one_of(all.vars(faceting)),
                screening, time_in_iso,
                contains("%")) %>%
    group_by(stringency, index_test_delay) %>%
    group_split %>%
    map(summarise_results, reduction = reduction) 
}



ribbon_plot <-
  function(x, 
           y_labels   = NULL, 
           colour_var = "stringency",
           by_type    = FALSE,
           custom_facets = NULL,
           bars  =TRUE
  ){
    
    
    if (is.null(custom_facets)){
      f_lhs <- c("waning", "index_test_delay", "delay_scaling")
      f_rhs <- c("yvar","stringency")
      
    } else {
      f_lhs <- all.vars(lhs(custom_facets))
      f_rhs <- all.vars(rhs(custom_facets))
    }
    
    if (!by_type){
      x <- filter(x, type == "all")
    } else {
      f_lhs <- c(f_lhs, "type") 
    }
    
    if (!any(f_rhs == "yvar") & length(y_labels) > 0){
      f_rhs <- f_rhs <- c("yvar", f_rhs) 
    }
    
    if (!all(f_lhs == ".")){
      f_lhs <- grep(pattern = ".", x = f_lhs, fixed = T, invert = T, value = T)
    }
    
    
    if (!is.null(y_labels)){
      x <- filter(x, yvar %in% names(y_labels))
      x <- mutate(x, yvar = factor(yvar, levels = names(y_labels), ordered = T))
    }
    
    # here we want to drop anything that's only got one value
    
    drop_unused_terms <- function(y,x){
      
      map_chr(.x = y, .f = function(v){
        if (length(unique(x[[v]])) > 1L | v == "."){
          v
        } else {NA_character_}
      }) %>% na.omit %>% c %>% unique
    }
    
    f_lhs <- drop_unused_terms(f_lhs, x)
    f_rhs <- drop_unused_terms(f_rhs, x)
    
    
    faceting_new <-
      as.formula(
        paste(
          paste(f_lhs,
                collapse = " + "),
          paste(f_rhs,
                collapse = " + "),
          sep = " ~ "
        )
      )
    
    
    x %<>% 
      test_labeller() %>% 
      mutate(type       = capitalize(type))
    
    xlims <- range(x$time_since_exp)
    
    colour_var_sym <- sym(colour_var)
    
    the_plot <- 
      ggplot(data = x, aes(x = time_since_exp,
                           y = M,
                           #color = !!colour_var,
                           fill  = !!colour_var_sym)) +
      facet_nested(nest_line = TRUE,
                   drop      = TRUE,
                   facets    = faceting_new,
                   labeller  = labeller(index_test_delay = index_test_labeller,
                                        delay_scaling    = 
                                          function(x){delay_scaling_labeller(x)},
                                        waning           = waning_labeller,
                                        #stringency       = test_labeller,
                                        type             = capitalize,
                                        yvar             = infectivity_labels,
                                        .multi_line      = TRUE)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            panel.border = element_rect(fill=NA),
            axis.ticks = element_line()) + 
      xlab(expression("Quarantine required until"~italic("t")~"days have passed since exposure")) +
      ylab("Median transmission potential averted") 
    
    if(bars == TRUE){
      the_plot <-  the_plot + 
        geom_linerange(aes(ymin = `2.5%`,
                           ymax = `97.5%`,
                           colour=!!colour_var_sym),size=1,alpha=0.3) +
        geom_linerange(aes(ymin = `25%`,
                           ymax = `75%`,
                           colour=!!colour_var_sym),size=1,alpha=0.5)
    } else {
      the_plot <- the_plot +
        geom_line(aes(y=`2.5%`,colour=!!colour_var_sym),linetype="dashed")+
        geom_line(aes(y=`97.5%`,colour=!!colour_var_sym),linetype="dashed")
    }
    
    the_plot <- the_plot +
      geom_point(aes(y = `50%`,
                     color = !!colour_var_sym),
                 pch="-",
                 size=5) +
      scale_x_continuous(minor_breaks = seq(xlims[1], xlims[2], by = 1),
                         breaks       = seq(xlims[1], xlims[2], by = 7))+
      scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))
    
    
    if (colour_var == "stringency"){
      the_plot <- the_plot +
        scale_color_manual(name = "Number of tests required before release",
                           values = covid_pal) +
        scale_fill_manual(name = "Number of tests required before release",
                          values = covid_pal)
    } else{ 
      the_plot <- the_plot + 
        scale_color_manual(name = "Type of infection",
                           values = covid_pal2) +
        scale_fill_manual(name = "Type of infection",
                          values = covid_pal2)
    } 
    
    the_plot
  } 

#' log scale
#'
#' Creates a function which returns ticks for a given data range. It uses some
#' code from scales::log_breaks, but in contrast to that function it not only
#' the exponentials of the base b, but log minor ticks (f*b^i, where f and i are 
#' integers), too.
#'
#' @param n Approximate number of ticks to produce
#' @param base Logarithm base
#'
#' @return
#'
#' A function which expects one parameter:
#'
#' * **x**: (numeric vector) The data for which to create a set of ticks.
#'
#' @export
logTicks <- function(n = 5, base = 10){
  # Divisors of the logarithm base. E.g. for base 10: 1, 2, 5, 10.
  divisors <- which((base / seq_len(base)) %% 1 == 0)
  mkTcks <- function(min, max, base, divisor){
    f <- seq(divisor, base, by = divisor)
    return(unique(c(base^min, as.vector(outer(f, base^(min:max), `*`)))))
  }
  
  function(x) {
    rng <- range(x, na.rm = TRUE)
    lrng <- log(rng, base = base)
    min <- floor(lrng[1])
    max <- ceiling(lrng[2])
    
    tck <- function(divisor){
      t <- mkTcks(min, max, base, divisor)
      t[t >= rng[1] & t <= rng[2]]
    }
    # For all possible divisors, produce a set of ticks and count how many ticks
    # result
    tcks <- lapply(divisors, function(d) tck(d))
    l <- vapply(tcks, length, numeric(1))
    
    # Take the set of ticks which is nearest to the desired number of ticks
    i <- which.min(abs(n - l))
    if(l[i] < 2){
      # The data range is too small to show more than 1 logarithm tick, fall
      # back to linear interpolation
      ticks <- pretty(x, n = n, min.n = 2)
    }else{
      ticks <- tcks[[i]]
    }
    return(ticks)
  }
}


# colours based on https://github.com/kgostic/traveller_screening/blob/7e7e84979b8c665392009d639dd9b3f2fc9b009a/nCov_manuscript_plots.R#L46

detected_levels <- c(#"Ceased being infectious before flight" = NA,
                     #"Develops symptoms before flight" = NA,
                     "Detected by pre-flight test"         = "#9ecae1", # detected pre-
                     "Detected by post-flight test"        = "#3182bd",
                     
                     "Develops symptoms on flight"         = "#ffffd4",
                     "Develops symptoms in quarantine"     = "#fee391",
                     "Develops symptoms in community"      = "#fec44f", 
                     "Not detected"                        = "#de2d26")      # detected post-)


infectious_levels <- c("Prevented from boarding"          = 0.2,
                       "Never infectious"                 = 0.33,
                       "No longer infectious"             = 0.5,
                       "Becomes infectious after release" = 0.67,
                       "Infectious due to low adherence"  = 0.8,
                       "Infectious when released"         = 1)

infectious_levels_edge <- c("Never infectious"                 = NA,
                            "No longer infectious"             = NA,
                            "Becomes infectious after release" = NA,
                            "Infectious when released"         = "black",
                            "Prevented from boarding"          = NA,
                            "Infectious due to low adherence"  = NA)

calculate_detected_infectious_status <- function(x, 
                                                 detected_levels, 
                                                 infectious_levels, 
                                                 drop_never_infectious = TRUE){
  
  x_ <- ungroup(x) %>%
    mutate(test_iso_end_t    = ifelse(is.infinite(test_iso_end_t), NA, test_iso_end_t),
           symp              = type == "symptomatic",
           symp_end_t        = ifelse(!symp, NA, symp_end_t),
           inf_dur           = pmax(inf_end,flight_dep_t) - pmax(inf_start,flight_dep_t),
           symp_post_arrival = symp & (onset_t >= flight_arr_t),
           symp_pre_arrival = symp & (onset_t <= flight_dep_t),
           symp_in_flight    = symp & (onset_t >= flight_dep_t) & (onset_t<=flight_arr_t),
           #symp_on_arrival   = symp_pre_arrival & !fly,
           symp_in_quar      = symp_post_arrival & (onset_t < quar_end_t), # quar doesn't cover daily testing
           symp_in_community = symp_post_arrival & (is.na(quar_end_t)), # does onset happen during daily testing?
           symp_after_quar   = symp & (onset_t >= quar_end_t),
           
           latest_release    = pmax(test_iso_end_t, quar_end_t, symp_end_t, na.rm = T),
           inf_after_release = inf_end >= latest_release
    ) %>%
    mutate(inf_after_release = replace_na(inf_after_release, TRUE)) %>% 
    mutate(detected =
             case_when(
               inf_dur==0           ~ "Ceased being infectious before flight",
              # symp_label           ~ "Symptomatic at departure", # ok
               pre_board_label      ~ "Detected by pre-flight test", # ok
               symp_pre_arrival     ~ "Develops symptoms before flight",
               symp_in_flight       ~ "Develops symptoms on flight",
               test_no != "None"    ~ "Detected by post-flight test", # daily test only? no! also contains quarantine with PCR/LFA testing
               symp_in_community    ~ "Develops symptoms in community",
               symp_in_quar & adhering_quar  ~ "Develops symptoms in quarantine",
               symp_in_quar & !adhering_quar  ~ "Develops symptoms in community",
               symp_after_quar      ~ "Develops symptoms in community",
               TRUE                 ~ "Not detected"
             ),
           infectious = 
             case_when(
              # symp_label                      ~ "Prevented from boarding", # ok
               pre_board_label                 ~ "Prevented from boarding", # ok
               is.infinite(inf_dur)            ~ "Never infectious",
               symp_in_quar   & !adhering_quar ~ "Infectious due to low adherence",
               symp_in_flight & !adhering_quar ~ "Infectious due to low adherence",
               grepl(x = detected, pattern = "^(Detected|Develops)") &
                 !adhering_symp|!adhering_test ~ "Infectious due to low adherence",
               inf_start >= latest_release     ~ "Becomes infectious after release",
               inf_after_release               ~ "Infectious when released",
               TRUE                            ~ "No longer infectious"
             )) %>%
    ungroup() 
    
    x_ %<>%
    mutate(detected    = factor(detected,   levels = names(detected_levels), ordered = T),
           infectious  = factor(infectious, levels = names(infectious_levels)))
  
  if (drop_never_infectious){
  x_ %<>% 
    filter(detected   != "Ceased being infectious before flight",
           infectious != "Never infectious")
  }
  
  x_
}

plot_outcomes <- function(x, v, 
                          xlab = "",
                          ylab = "Proportion of infectious travellers",
                          faceting = ~ stringency + pre_board_assay){
  
  # need to recalculate the AUC curves so that we get the right denominator for 
  # source of infectivity

  x_to_plot <- x %>%
    filter(!is.na({{ v }})) %>%
    mutate(infectious = ifelse(grepl(x = infectious,
                                     pattern = "(I|i)nfectious.*(release|adherence)"),
                               "Risk of transmission in community",
                               "No risk of transmission in community")) %>%
    mutate(detected = factor(detected, levels = names(detected_levels))) %>%
    arrange(detected) %>%
    mutate(group = factor(interaction(detected, infectious, sep = " "))) %>%
    arrange(infectious, group) %>%
    mutate(group = fct_inorder(group)) 

  ggplot(data = x_to_plot,
         aes(x = factor({{ v }}))) +
    geom_bar(aes(fill   = detected,
                 group  = group,
                 alpha  = infectious
                 ),
             position = position_fill()) +
    geom_bar(aes(color  = infectious),
             fill = NA,
             #color  = infectious),
             position = position_fill()) +
    facet_nested(faceting,
                 nest_line = TRUE,
                 labeller = labeller(
                   stringency     = label_value,
                   sens_LFA       = sens_scaling_labels,
                   type           = type_labels,
                   pre_board_test = c(`TRUE` = "Pre-flight test",
                                      `FALSE` = "No pre-flight test"),
                   assay = function(x){paste(x, "post-arrival")},
                   adherence_quar = function(x){label_template(x, "Quar. Ad:")},
                   adherence_symp = function(x){label_template(x, "Symp. Ad:")},
                   adherence_test = function(x){label_template(x, "Pos.-test Ad:")},
                   adherence_iso  = function(x){label_template(x, "Iso. Ad:")})) + 
    scale_fill_manual(values = detected_levels,
                      breaks = names(detected_levels),
                      drop   = FALSE,
                      name   = "Method of detection") +
    scale_alpha_manual(values = c(0.5, 0.8),
                       name   = "Infectivity status",guide=F) +
    scale_colour_manual(values = c(NA,"black"),
                        name   = "Infectivity status")+
    plotting_theme +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(size = 6)) +
    labs(x = xlab, y = ylab) + 
    guides(fill  = guide_legend(nrow = 2, byrow=TRUE),
           # alpha = guide_legend(nrow = 1, byrow = TRUE,
           #                      override.aes = list(fill = c(alpha("grey40", 0.2),
           #                                                  "grey40"))
           #                      ),
           color = guide_legend(nrow = 1, byrow = TRUE)) +
    #alpha = guide_legend(override.aes = list(fill = "black"))) +
    scale_y_continuous(breaks = breaks_width(0.25),#round(x = seq(0,1,length.out = 7), 2),
                       labels = function(x){round(x,2)}) 
}




make_results_plot <- function(x, faceting = ~stringency + pre_board_assay, ...){
  # ... are arguments to save_plot
  
  daily <- plot_outcomes(x = x,
                         v = n_tests,
                         xlab = "Number of days of testing",
                         faceting = faceting)
  
  quar  <- plot_outcomes(x = x,
                         v = quar_dur,
                         xlab = "Duration of quarantine (days)",
                         faceting = faceting) 
  
  col_facets <- formula.tools::rhs.vars(faceting)
  row_facets <- formula.tools::lhs.vars(faceting)
  
  widths <- x %>% 
    select(quar_dur, n_tests, all_of(col_facets)) %>% # select all variables that contribute a bar
    distinct %>%
    gather(key, value, n_tests, quar_dur) %>%
    mutate(key = factor(key, levels = c("quar_dur", "n_tests"))) %>%
    na.omit %>%
    count(key) %>% 
    pull(n) %>%
    {. + pmax(2, length(row_facets))}  # need some space for axis labels and faceting
  
  
  # 40, 12 is verrrry rigid. should be based on how many unique values of x we see!
  # outcomes_plot <- quar + daily + plot_layout(nrow = 1, 
  #                                             widths = widths,   
  #               guides = 'collect') &
  #   theme(legend.position = "bottom",
  #         legend.box = "vertical")
  
  theme_bit <- theme(legend.position="none",
          plot.margin = margin(6, 4, 6, 4))
  
  leg <- cowplot::get_legend(quar + theme(legend.box.margin = margin(12, 6, 4, 6)))
  p   <- cowplot::plot_grid(quar  + theme_bit , 
                            daily + theme_bit,
                            rel_widths = widths, labels = c("A", "B"),
                            nrow = 1)
  
  outcomes_plot <- cowplot::plot_grid(p, leg, ncol = 1, rel_heights = c(1, 0.25))
  
  
  
  
  save_plot(plot = outcomes_plot, ...)
  
}



axis_label_lt <- function(x){
  
  s <- scales::scientific_format()(x) %>%
    sub(pattern = "e\\+",
        x = .,
        replacement = "e")
  
  index <- min(which(!is.na(x)))
  
  s[index] <- paste0("''<=", s[index])
  
  s <- gsub(pattern = "1e", replacement = "10^", x = s)
  
  s <- gsub(pattern = "e",
            replacement = " %*% 10^",
            x = s)
  parse(text = s)
  
}

axis_label_lt_percent <- function(x, ...){
  
  s <- prettyNum(sprintf("%g", 100*x), big.mark = ",")
  #s %<>% gsub(pattern = "(\\d+)\\%", replacement = )
  index <- min(which(!is.na(x)))
  s[index] <- paste0("$\\leq ", s[index])
  s <- paste0(s, "\\%$")
  
  latex2exp::TeX(s)
  
}



trim_zeros <- function(x){
  #browser()
  comma(x) %>% 
    str_remove(string = ., pattern =  "0+$") %>%
    str_remove(string = ., pattern = "\\.$")
}

label_number_si2 <- function(x,  lt = FALSE, ...){
  
  accuracy <- na.omit(x) %>% min %>% log10 %>% floor %>% {10^.}
  
  s <- label_number_si(accuracy = accuracy, drop0trailing = TRUE, ...)(x)
  
  if (lt){
    index <- min(which(!is.na(x)))
    
    s[index] <- paste0("\\leq ", s[index])
    
    s <- paste0("$", s, "$")
    
    return(latex2exp::TeX(s))
  } else {
    return(latex2exp::TeX(s))
  }
  
}


traj_pal <- RColorBrewer::brewer.pal(name="Dark2",n=4)[2:4]
names(traj_pal) <- c("LFT","PCR","Culture")

post_test_pal <- c("No testing" = "#808080",
                   traj_pal)



xrange <- data.frame(xmin = 0.01, xmax = 1e6) %>% gather(key, value)
# 
# fill_3 <- list(`Less than 1%`      = mutate(xrange, ymin = 0,    ymax = 0.01),
#                `1%-10%`            = mutate(xrange, ymin = 0.01, ymax = 0.1 ),
#                `10%-100%`          = mutate(xrange, ymin = 0.1,  ymax = 1   ),
#                `Greater than 100%` = mutate(xrange, ymin = 1  ,  ymax = Inf )) %>%
#   bind_rows(.id = "band") %>%
#   mutate_at(.vars = vars(ymin, ymax), .funs = ~{. * value}) %>%
#   mutate(band = fct_inorder(band)) %>%
#   crossing(variant = c(TRUE, FALSE))

fill_4 <- list(`Less than 1%`      = c(L = 0,    U = 0.01),
               `1%-10%`            = c(L = 0.01, U = 0.1),
               `10%-100%`          = c(L = 0.1,  U = 1),
               `Greater than 100%` = c(L = 1,    U = Inf)) %>%
  map_df(~mutate(xrange, 
                 ymin = .x['L'] ,
                 ymax = .x['U'] ),
         .id = 'band') %>%
  mutate(band = fct_inorder(band))

fill_3 <- fill_4 %>% inner_join(xrange) %>%
  #mutate_at(.vars = vars(ymin, ymax), .funs = ~{. / value}) %>%
  crossing(variant = c(TRUE, FALSE))

shapes <- c(
  "Pre-flight LFT + daily testing (5)"     = 23,
  "Pre-flight LFT + quarantine (5) + LFT"  = 25,
  "Pre-flight LFT + quarantine (10) + LFT" = 24, 
  "Pre-flight LFT"                         = 22,
  "No intervention"                        =  4)

country_list <- c("USA", "ISR", "LUX", "MEX", "RUS", "SGP")
