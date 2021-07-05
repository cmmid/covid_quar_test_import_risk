pacman::p_load(readr,brms,broom,modelr,patchwork,tidyverse,pROC,epiR,pdftools,fitdistrplus,conflicted,ragg,here)
conflicted::conflict_prefer(name = "select",winner="dplyr")
conflicted::conflict_prefer(name = "filter",winner="dplyr")

##### KCL ANALYSIS ----


pickering <- readxl::read_xlsx(here("data","pickering_dat.xlsx")) %>%
  select(-c(`Viral Growth`,...7,...8)) %>%
  rename("Culture"=...6) %>%
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)ifelse(x=="ND",NA,x)) %>%
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)case_when(x%in%c(0.5,1,2)~1,
                                         is.na(x)~NA_real_,
                                         TRUE~0)) %>%
  mutate(id=row_number()) %>%
  rename(ct=`Ct N1`)


pickering %>% 
  pivot_longer(cols = c(`SureScreen F`,Innova,Encode,Culture)) %>% 
  ggplot()+geom_jitter(aes(x=ct,y=name,colour=factor(value)),size=5,alpha=0.5,width = 0,height =0.2)


pickering %>% 
  rename("Innova LFT"=Innova,"Viral culture"=Culture) %>% 
  pivot_longer(cols = c(`SureScreen F`,`Innova LFT`,Encode,`Viral culture`)) %>% 
  filter(name%in%c("Viral culture","Innova LFT")) %>% 
 # bind_rows(innova_liv_sim %>% select(ct,AG) %>% rename("value"=AG) %>% mutate(name="Liverpool")) %>% 
  drop_na(value) %>% 
  ggplot(aes(x=ct,y=value,colour=name))+
  geom_jitter(size=2,alpha=0.5,width = 0,height =0.05)+
  geom_smooth(method = "glm",se=F,method.args=list(family="binomial"),size=1)+ 
  # geom_line(data=infectivity %>% select(ct,value) %>% 
  #             mutate(name="Lee et al. infectivity"),
  #           aes(x=ct,y=value,colour=name),size=1)+
  geom_line(data=tribble(~ct,~value,
                         10,1,
                         35,1,
                         35,0,
                         40,0) %>% 
              mutate(name="PCR"),
            aes(x=ct,y=value,colour=name),
            size=1) +
  scale_color_brewer(type="qual",palette = "Set1",direction=-1)+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA))+
  scale_x_continuous(breaks=scales::breaks_width(5))+
  #scale_x_continuous(breaks=seq(0,9),labels=c("0\n37.1","1\n33.8","2\n30.5", "3\n27.2", "4\n23.9", "5\n20.6", "6\n17.35", "7\n14.1", "8\n10.8",  "9\n7.5" ))+
  labs(x="Cycle Threshold",
       y="Probability of detection/ infectivity",
       colour="")

save_plot(dpi = 400, 
          device = "png",
          prefix = "Kings_detection",
          base = "plot", 
          width = 210, 
          height = 120)

innova_mod  <- glm(Innova~ct,data=pickering,family="binomial") 
culture_mod <- glm(Culture~ct,data=pickering,family="binomial") 


pickering %>% mutate(name="Innova KCL") %>%
  ggplot(aes(x=ct,y=Innova,colour=name))+
  #geom_point(size=1,alpha=0.5,width = 0,height =0.2)+
  geom_smooth(method = "glm",se=F,method.args=list(family="binomial"),size=1)+ 
  geom_line(data=tribble(~ct,~Innova,
                                10,1,
                                35,1,
                                35,0,
                                40,0) %>% mutate(name="PCR"),
            aes(x=ct,y=Innova,colour=name),
            size=1) +
  scale_color_manual(values=c("#D95F02","#7570B3","#E7298A","#66A61E"))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA))+
  scale_x_continuous(breaks=scales::breaks_width(5))+
  labs(x="Cycle Threshold",
       y="Probability of detection",
       colour="")
