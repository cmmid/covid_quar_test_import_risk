## curves for paper

# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")
source("lfa_test_data.R")


set.seed(2020)

traj <- readRDS(file.info(list.files("results/", full.names = T)) %>% 
                  as.data.frame() %>% 
                  rownames_to_column()%>% 
                  filter(str_detect(rowname,"traj_")) %>% 
                  slice_max(mtime) %>% 
                  pull(rowname))

trajectories_to_plot <- traj %>%
  #filter(type=="symptomatic") %>% 
  group_by(type) %>% 
  sample_n(4) %>% 
  #ungroup() %>% 
  mutate(pred = map(.x = m, 
                    ~data.frame(x = seq(0, 25, length.out = 101)) %>%
                      mutate(ct = .x(x)))) %>%
  unnest(pred) %>%
  mutate(LFT = stats::predict(innova_mod,
            newdata =
              data.frame(ct = ct),
            type = "response"
  ),
  culture=stats::predict(culture_mod,
            newdata =
              data.frame(ct = ct),
            type = "response"
  ),
  PCR=ifelse(ct<35,1,0)) %>% 
  arrange(type)

traj_a <- trajectories_to_plot %>% 
  #pivot_longer(cols=c(ct,culture,LFT)) %>%  
  ggplot()+
  geom_line(aes(x=x,
                y=ct,
                group=interaction(sim,idx),
                colour="Cycle threshold",
                size=(factor(type))))+
  scale_color_manual(values="black",name="",guide=F)+
  scale_x_continuous(name="Days since exposure",breaks = breaks_width(7))+
  scale_size_manual(name="",labels=capitalize, values = c("symptomatic" = 1, 'asymptomatic' = 0.5))+
  scale_y_reverse(# Features of the first axis
    limits=c(39.99,NA),
    name = "Cycle\nthreshold")+
  plotting_theme+
  theme( strip.background = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text.x = element_blank(),
         #axis.title.x = element_blank()
         )+
  facet_wrap(~fct_rev(type)+interaction(sim,idx),nrow = 1)

traj_b <- trajectories_to_plot %>% 
  #pivot_longer(cols=c(ct,culture,LFT)) %>%  
  ggplot()+
  geom_line(aes(x=x,
                y=culture,
                group=interaction(sim,idx),
                colour="Culture",
                size=fct_rev(type)))+
  geom_line(aes(x=x,
                y=LFT,
                group=interaction(sim,idx),
                colour="LFT",
                size=fct_rev(type)))+
  geom_line(aes(x=x,
                y=PCR,
                group=interaction(sim,idx),
                colour="PCR",
                size=fct_rev(type)))+
  scale_color_manual(values=traj_pal,name="")+
  scale_size_manual(name="",labels=capitalize, guide = F,
                    values = c("symptomatic" = 1, 'asymptomatic' = 0.5))+
  scale_x_continuous(name="Days since exposure",breaks = breaks_width(7))+
  scale_y_continuous(# Features of the first axis
    name = "Probability\nof detection")+
  plotting_theme+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank()) +
  facet_wrap(~fct_rev(type)+interaction(sim,idx),nrow = 1) +
  guides(color = guide_legend(override.aes = list(size = 2) ) )

pickering_plot <- pickering %>% 
  rename("LFT"=Innova,"Culture"=Culture) %>% 
  pivot_longer(cols = c(`SureScreen F`,`LFT`,Encode,`Culture`)) %>% 
  filter(name%in%c("Culture","LFT")) %>% 
  # bind_rows(innova_liv_sim %>% select(ct,AG) %>% rename("value"=AG) %>% mutate(name="Liverpool")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x=ct,y=value,colour=name))+
  #geom_jitter(alpha=0.25,width = 0, height =0.05)+
  geom_point(alpha = 0.25, pch = "|", size = 4) +
  geom_smooth(method = "glm",se=F,method.args=list(family="binomial"),size=0.5)+ 
  # geom_line(data=infectivity %>% select(ct,value) %>% 
  #             mutate(name="Lee et al. infectivity"),
  #           aes(x=ct,y=value,colour=name),size=1)+
  geom_line(data=tribble(~ct,~value,
                         10,1,
                         35 ,1,
                         35.0001,0,
                         40,0) %>% 
              mutate(name="PCR"),
            aes(x=ct,y=value,colour=name),
            size=0.5) +
  scale_color_manual(values=traj_pal,name="",guide=F)+
  plotting_theme + 
  theme(panel.border = element_rect(fill=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_x_reverse(
    #breaks=scales::breaks_width(5)
    )+
  #scale_x_continuous(breaks=seq(0,9),labels=c("0\n37.1","1\n33.8","2\n30.5", "3\n27.2", "4\n23.9", "5\n20.6", "6\n17.35", "7\n14.1", "8\n10.8",  "9\n7.5" ))+
  labs(x="Cycle threshold",
       y="Probability\nof detection",
       colour="")+
  facet_wrap(~name)


trajectories_plot <- traj_a/pickering_plot/traj_b+ 
  plot_layout(ncol=1,heights = c(1,1.5,1),guides="collect")&
  plot_annotation(tag_levels = "A")&
  theme(legend.position = "bottom",strip.placement = NULL)

save_plot(plot = trajectories_plot, dpi=400,
          device="png",prefix = "trajectories", height = 150,width=210)

trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=x,
                y=ct,
                group=interaction(sim,idx)),
            colour=traj_pal[1],show.legend = F)+
  # geom_line(aes(x=x,
  #               y=culture*40,
  #               group=interaction(sim,idx),
  #               colour="Culture"))+
  # geom_line(aes(x=x,
  #               y=LFT*40,
  #               group=interaction(sim,idx),
  #               colour="LFT"))+
  scale_colour_viridis_d(name="",guide=F)+
  scale_x_continuous(name="Days since exposure")+
  scale_y_continuous(name="CT")+
  theme_half_open()

save_plot(dpi=400,
          device="png",prefix = "trajectories_ct", height = 100,width=100)

trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=x,
                y=culture,
                group=interaction(sim,idx)),
                colour=traj_pal[2],show.legend = F)+
  # geom_line(aes(x=x,
  #               y=LFT*40,
  #               group=interaction(sim,idx),
  #               colour="LFT"))+
  scale_x_continuous(name="Days since exposure")+
  scale_y_continuous(name="Culture probability",limits=c(0,1))+
  theme_half_open()

save_plot(dpi=400,
          device="png",prefix = "trajectories_culture", height = 100,width=100)

trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=x,
                y=LFT,
                group=interaction(sim,idx)),
            colour=traj_pal[3],show.legend = F)+
  # geom_line(aes(x=x,
  #               y=LFT*40,
  #               group=interaction(sim,idx),
  #               colour="LFT"))+
  scale_x_continuous(name="Days since exposure")+
  scale_y_continuous(name="LFT probability",limits=c(0,1))+
  theme_half_open()

save_plot(dpi=400,
          device="png",prefix = "trajectories_lft", height = 100,width=100)

trajectories_to_plot %>% 
  ggplot()+
  geom_line(aes(x=x,
                y=PCR,
                group=interaction(sim,idx)),
            colour=traj_pal[4],show.legend = F)+
  # geom_line(aes(x=x,
  #               y=LFT*40,
  #               group=interaction(sim,idx),
  #               colour="LFT"))+
  scale_x_continuous(name="Days since exposure")+
  scale_y_continuous(name="PCR probability",limits=c(0,1))+
  theme_half_open()

save_plot(dpi=400,
          device="png",prefix = "trajectories_pcr", height = 100,width=100)


# lft_plot <- trajectories_to_plot %>% 
#   ggplot(data = ., aes(x = x, y = y)) +
#   # annotate(geom = "rect",
#   #          ymax = Inf, xmin = -Inf, xmax = Inf, ymin = 30,
#   #          color = NA,
#   #          fill = "black", alpha = 0.1) +
#   geom_line(aes(group = idx,
#                 color = LFT,
#                 ),
#             alpha = 0.1) +
#   theme_bw() +
#   facet_wrap(~value,labeller = labeller(value=capitalize),ncol=1)+
#   #geom_hline(aes(yintercept=35,linetype="PCR detection threshold"))+
#   #geom_hline(aes(yintercept=30,linetype="Infectivity threshold"))+
#   scale_linetype_manual(name="",values=c("dotted","dashed"))+
#   # geom_hline(data = mask,
#   #            aes(yintercept = max),
#   #            lty = 2) +
#    scale_color_viridis_c(option="magma",begin=0.2,end=0.8,
#                       name   = "Probability of detection by LFT",guide=guide_colorsteps(ticks=T,barwidth=unit(5,"cm"))) +
#   ylab("Ct value") +
#   xlab("Time since exposure (days)") +
#   plotting_theme +
#   theme(panel.spacing = unit(1,"cm")) +
#   scale_y_reverse()+coord_cartesian(expand = F)
# 
# 
# culture_plot <- trajectories_to_plot %>% 
#   ggplot(data = ., aes(x = x, y = y)) +
#   geom_line(aes(group = idx,
#                 color = culture,
#   ),
#   alpha = 0.1) +
#   theme_bw() +
#   facet_wrap(~value,labeller = labeller(value=capitalize),ncol=1)+
#   scale_linetype_manual(name="",values=c("dotted","dashed"))+
#   scale_color_viridis_c(option="viridis",begin=0.2,end=0.8,
#                         name   = "Probability of viral culture",guide=guide_colorsteps(ticks=T,barwidth=unit(5,"cm"))) +
#   ylab("Ct value") +
#   xlab("Time since exposure (days)") +
#   plotting_theme +
#   theme(panel.spacing = unit(1,"cm")) +
#   scale_y_reverse()+coord_cartesian(expand = F)
# 
# lft_plot+culture_plot
# 
# save_plot(dpi=400,
#           device="png",prefix = "trajectories", height = 200,width=210)
# 
# 
# trajectories_to_plot %>%
#   ggplot(data = ., aes(x = x, y = culture)) +
#   geom_line(aes(group = idx,
#                 #color = culture,
#   ),
#   alpha = 0.1) +
#   stat_summary(fun="mean",geom="line",colour="red",linetype="dashed")+
#   stat_summary(fun="median",geom="line",colour="red",linetype="solid")+
#   theme_bw() +
#   facet_wrap(~value,labeller = labeller(value=capitalize),ncol=1)+
#   scale_linetype_manual(name="",values=c("dotted","dashed"))+
#   scale_color_viridis_c(option="viridis",begin=0.2,end=0.8,
#                         name   = "Probability of viral culture",guide=guide_colorsteps(ticks=T,barwidth=unit(5,"cm"))) +
#   ylab("Probability of culture") +
#   xlab("Time since exposure (days)") +
#   plotting_theme +
#   theme(panel.spacing = unit(1,"cm")) +
#   #scale_y_reverse()+
#   coord_cartesian(expand = F)
# 
