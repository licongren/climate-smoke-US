library(openxlsx)
library(ggplot2)
library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(MetBrewer)
library(dplyr)
library(purrr)
library(data.table)
library(cowplot)
library(zoo)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/health"

setwd(figure_path)
states_map <- map_data("state")

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"

#-------------------------------------------------------------------------------
# Plot the smoke-mortality relationship
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------


##----- linear 
coef <- readRDS(paste0(result_path, "/coef_linear.rds")) %>% filter(age_group!="under_5")

coef_plot <- coef %>% group_by(age_group) %>%
  summarise(median=quantile(coef,0.5),
            p025=quantile(coef,0.025),
            p975=quantile(coef,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=age_group)) +
  geom_point(aes(x=age_group, y=median), position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=p025,
                    ymax=p975), position=position_dodge(width=0.2), width=0.01) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(x="", y="mortalities per 100k\nper 1ug smoke PM2.5") +
  geom_hline(aes(yintercept=0), linetype="dashed")
ggsave("linear_annual_boot.pdf", width=5, height=3.5)

##----- linear bins
coef <- readRDS(paste0(result_path, "/coef_linear_bins.rds")) %>% filter(age_group!="under_5")

coef_plot <- coef %>% filter(panel=="year", model=="linear bins (coarse)") %>%
  group_by(age_group, bins) %>%
  summarise(median=quantile(coef,0.5),
            p025=quantile(coef,0.025),
            p975=quantile(coef,0.975)) %>% ungroup()
  
coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=bins, y=median)) +
  geom_point()+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975),
                width=0.001) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=14))
ggsave("linear_year_coarse_bins_boot.pdf", width=5.5, height=6)

### fine bins
coef_plot <- coef %>% filter(panel=="year", model=="linear bins") %>%
  group_by(age_group, bins) %>%
  summarise(median=quantile(coef,0.5),
            p025=quantile(coef,0.025),
            p975=quantile(coef,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=bins, y=median)) +
  geom_point()+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975),
                width=0.001) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=14))
ggsave("linear_year_bins_boot.pdf", width=7, height=6)
 
##----- quadratic and cubic: create response functions
smoke_conc_vec <- c(seq(0,2.99, by=0.01), seq(3,10,by=0.1))
nboot <- 500

coef <- readRDS(paste0(result_path, "/coef_quadratic.rds")) %>% filter(age_group!="under_5")
response <- NULL
for (aa in unique(coef$age_group)){
  for (bbb in 1:nboot){
    print(bbb)
    coef_tmp <- filter(coef, bootid==bbb, age_group==aa)
  response_tmp <- data.frame(smoke=smoke_conc_vec) %>%
                   mutate(death_rate=as.numeric(coef_tmp[1])*smoke+as.numeric(coef_tmp[2])*smoke^2,
                          bootid=bbb, age_group=aa, model="quadratic")
  response <- bind_rows(response_tmp, response)
  }
}

coef <- readRDS(paste0(result_path, "/coef_cubic.rds")) %>% filter(age_group!="under_5")
for (aa in unique(coef$age_group)){
  for (bbb in 1:nboot){
    print(bbb)
    coef_tmp <- filter(coef, bootid==bbb, age_group==aa)
    response_tmp <- data.frame(smoke=smoke_conc_vec) %>%
      mutate(death_rate=as.numeric(coef_tmp[1])*smoke+as.numeric(coef_tmp[2])*smoke^2+as.numeric(coef_tmp[3])*smoke^3,
             bootid=bbb, age_group=aa, model="cubic")
    response <- bind_rows(response_tmp, response)
  }
}
saveRDS(response, paste0(result_path, "/quadratic_response_curve.rds"))

response_summ <- response %>% group_by(age_group, smoke, model) %>%
summarise(median=median(death_rate),
          p025=quantile(death_rate, 0.025),
          p975=quantile(death_rate, 0.975)) %>%
  ungroup()

ggplot(response_summ %>% filter(smoke<10)) +
  geom_line(aes(x=smoke, y=median, colour=model)) +
  geom_ribbon(aes(x=smoke, ymin=p025, ymax=p975, colour=model, fill=model), alpha=0.5) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=18)) +
  labs(x="", y="mortalities per 100k\nper 1ug smoke PM2.5") +
  scale_colour_manual(values=c("darkorange", "purple4")) +
  scale_fill_manual(values=c("darkorange", "purple4"))
ggsave("year_quadratic_cubic.png", width=5, height=6)



######################
##----- poisson 
#######################
coef <- readRDS(paste0(result_path, "/coef_poisson.rds")) %>% filter(age_group!="under_5")

coef_plot <- coef %>% group_by(age_group) %>%
  summarise(median=quantile(coef,0.5),
            p025=quantile(coef,0.025),
            p975=quantile(coef,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=age_group)) +
  geom_point(aes(x=age_group, y=median), position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=p025,
                    ymax=p975), position=position_dodge(width=0.2), width=0.01) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x="", y="mortalities rate change\nper 1ug smoke PM2.5") +
  geom_hline(aes(yintercept=0), linetype="dashed")
ggsave("poisson_annual_boot.pdf", width=5, height=3.5)

##----- poisson bins
coef <- readRDS(paste0(result_path, "/coef_poisson_bins.rds")) %>% filter(age_group!="under_5")

coef_plot <- coef %>% filter(panel=="year", model=="poisson bins (coarse)") %>%
  group_by(age_group, bins) %>%
  summarise(median=quantile(exp(coef)-1,0.5),
            p025=quantile(exp(coef)-1,0.025),
            p975=quantile(exp(coef)-1,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=bins, y=median)) +
  geom_point()+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975),
                width=0.001) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=14)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x="", y="mortalities rate change\nper 1ug smoke PM2.5") 
ggsave("poisson_year_coarse_bins_boot.pdf", width=5.5, height=6)

### fine bins
coef_plot <- coef %>% filter(panel=="year", model=="poisson bins") %>%
  group_by(age_group, bins) %>%
  summarise(median=quantile(exp(coef)-1,0.5),
            p025=quantile(exp(coef)-1,0.025),
            p975=quantile(exp(coef)-1,0.975)) %>% ungroup()

coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))

ggplot(coef_plot, aes(x=bins, y=median)) +
  geom_point()+
  geom_errorbar(aes(ymin=p025,
                    ymax=p975),
                width=0.001) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=14)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x="", y="mortalities rate change\nper 1ug smoke PM2.5") 
ggsave("poisson_year_bins_boot.pdf", width=7, height=6)



