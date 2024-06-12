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
# analysis of county smoke-mortality robustness
# Written by Minghao
# Last edited Feb 2024
#-------------------------------------------------------------------------------
mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))
age_group <- "all_ages"

annual_fe <- "fips + state_name^year"
month_fe <- "fips^month + state_name^year + fips + year^month"
nboot <- 500

#-------------------------------------------------------------------------
#------ Running poisson models on annual outcomes (including 2020) -------
#-------------------------------------------------------------------------
#### poisson model on annual conc bins
cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

### annual bins
fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:10,1]) %>% 
      mutate(age_group=aa,  model="include 2020", 
             bins=levels(annual_bin$smoke_bins)[2:11], bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

#saveRDS(coef_combined, paste0(result_path, "/coef_poisson_bins_include_2020.rds"))

#-------------------------------------------------------------------------
#------ Running poisson models on annual binned monthly conc -------
#-------------------------------------------------------------------------
#### poisson model on annual conc bins
cutoff <- c(-0.01,0,0.25,0.5,0.75,1,2,3,4,5,10,20,150)
month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          month_fe))

#coef_combined <- NULL
for (aa in age_group){
  df <- filter(month_bin, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:11,1]) %>% 
      mutate(age_group=aa, model="month bins", 
             bins=levels(month_bin$smoke_bins)[2:12], bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

#-------------------------------------------------------------------------
#-----------------        Running annual bin of # months  ---------------
#-------------------------------------------------------------------------
cutoff <- c(-0.01,0,0.25,0.5,0.75,1,2,3,4,5,10,20,150)
month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

annual_month_bin <- month_bin %>% group_by(fips, year, smoke_bins, age_group) %>%
  summarise(n=length(smoke_bins)) %>% ungroup() %>%
  pivot_wider(id_cols = fips:age_group, names_from = smoke_bins, values_from = n)
colnames(annual_month_bin)[4:15] <- c("smoke_0","smoke_0_025","smoke_05_075", 
                                      "smoke_025_05","smoke_075_1", "smoke_1_2", 
                                      "smoke_4_5","smoke_5_10", "smoke_3_4",
                                      "smoke_2_3", "smoke_10_20", "smoke_20_plus")

annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  left_join(annual_month_bin) %>%
  mutate(across(c("smoke_0","smoke_0_025","smoke_05_075", 
                  "smoke_025_05","smoke_075_1", "smoke_1_2", 
                  "smoke_4_5","smoke_5_10", "smoke_3_4",
                  "smoke_2_3", "smoke_10_20", "smoke_20_plus"), ~replace_na(.x, 0)))

fmla <- as.formula(paste0("rate ~ smoke_0_025 + smoke_025_05 + smoke_05_075 +
                                  smoke_075_1 + smoke_1_2 +  smoke_2_3 + smoke_3_4 +
                                  smoke_4_5 + smoke_5_10 + smoke_10_20 + smoke_20_plus + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

#coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:11,1]) %>% 
      mutate(age_group=aa, model="bins in monthly conc", 
             bins=levels(month_bin$smoke_bins)[2:12], bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

#-------------------------------------------------------------------------
#-----------------        Running annual bin of # days  ---------------
#-------------------------------------------------------------------------
##create daily county-level smoke
cutoff <- c(-0.01,0.5,1,2,3,4,5,10,20,50, 2000)
daily_bin <- readRDS(paste0(data_path, "/smokePM2pt5_predictions_daily_county_20060101-20201231.rds"))  %>%
  mutate(year=as.numeric(substr(date,1,4))) %>%
  mutate(fips=as.numeric(GEOID)) %>%
  mutate(smoke_bins=cut(smokePM_pred, cutoff)) %>% ungroup()

annual_daily_bin <- daily_bin %>% group_by(fips, year, smoke_bins) %>%
  summarise(n=length(smoke_bins)) %>% ungroup() %>%
  pivot_wider(id_cols = fips:year, names_from = smoke_bins, values_from = n)
colnames(annual_daily_bin)[3:12] <- c("smoke_0","smoke_05_1","smoke_1_2", 
                                      "smoke_2_3", "smoke_3_4", "smoke_4_5",
                                      "smoke_5_10", "smoke_10_20", "smoke_20_50", "smoke_50_plus") 
## note the # of days in the lowest smoke bins is not correct  (but we did not fix since it will de dropped anyway)

annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  left_join(annual_daily_bin) %>%
  mutate(across(c("smoke_0","smoke_05_1","smoke_1_2", 
                  "smoke_2_3", "smoke_3_4", "smoke_4_5",
                  "smoke_5_10", "smoke_10_20", "smoke_20_50", "smoke_50_plus") , ~replace_na(.x, 0)))

fmla <- as.formula(paste0("rate ~ smoke_05_1 +  smoke_1_2 +  smoke_2_3 + smoke_3_4 +
                                  smoke_4_5 + smoke_5_10 + smoke_10_20 + smoke_20_50 + smoke_50_plus + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

#coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:9,1]) %>% 
      mutate(age_group=aa, model="bins in daily conc", 
             bins=levels(daily_bin$smoke_bins)[2:10], bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

saveRDS(coef_combined, paste0(result_path, "/coef_poisson_bins_robust.rds"))


###------------------------
## Plot the coefficient
###-----------------------
coef_robust <- readRDS(paste0(result_path, "/coef_poisson_bins_robust.rds"))

coef_summ <- coef_robust %>% group_by(age_group, model, bins) %>%
  summarise(mean=mean(exp(coef)-1, na.rm=T),
            p025=quantile(exp(coef)-1,0.025, na.rm=T),
            p975=quantile(exp(coef)-1, 0.975), na.rm=T) %>% ungroup()
coef_summ$bins <- factor(coef_summ$bins,
                         levels=c("(0,0.25]", "(0.1,0.25]", "(0.25,0.5]", "(0.5,0.75]","(0.5,1]",
                                  "(0.75,1]", "(1,2]", "(2,3]",  "(3,4]", "(4,5]", 
                                  "(5,6]", "(6,10]", "(5,10]","(10,20]", "(20,50]", "(20,150]","(50,2e+03]"))  
coef_summ$model <- factor(coef_summ$model,
                           levels=c("include 2020", "month bins", "bins in daily conc", "bins in monthly conc"))
  
ggplot(coef_summ, aes(x=bins, colour=model)) +
  geom_point(aes(x=bins, y=mean), position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=p025,
                    ymax=p975), position=position_dodge(width=0.2), width=0.01) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(x="", y="mortalities rate changes ") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  facet_wrap(~model, scales = "free", ncol=1)
ggsave("poisson_bins_robust.png", width=11, height=9.5)

###------------------------------------------------------
## Calculate the moralities using robust poisson bin CRFs
###------------------------------------------------------
health_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/health"
smoke_obs <- readRDS(paste0(data_path, "/health/US_smoke_county_2011_2020.rds"))  %>%
  filter(age_group == "all_ages")
coef_robust <- readRDS(paste0(result_path, "/coef_poisson_bins_robust.rds"))
nboot <- max(coef_robust$bootid)

##------- Main estimates
main <- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins", scenario=="Observation") %>%
  mutate(spec="main")

alternative_bins<- readRDS(paste0(health_path, "/summ_smoke_death_smokePM_CRF_mediangcm.rds")) %>%
  filter(age_group=="all_ages", CRF=="poisson bins (coarse)", scenario=="Observation") %>%
  mutate(spec="alternative bins")
#--- include 2020
cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
include_2020 <- smoke_obs %>% 
  mutate(smokePM=replace(smokePM, smokePM>10, 10)) %>%
  mutate(bins=cut(smokePM, cutoff)) %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(smoke_obs))) %>%
  left_join(filter(coef_robust, model=="include 2020")) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = (exp(coef)-1)*death_rate_avg*pop_2019) %>%
  group_by(year, bootid) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  group_by(bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% mutate(spec="include 2020")

#--- month bins
cutoff <- c(-0.01,0,0.25,0.5,0.75,1,2,3,4,5,10,20,150)
month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(bins=cut(smokePM_mean, cutoff)) %>%
  ungroup() %>%
  group_by(fips, year, bins) %>%
  summarise(nbins=length(bins)) %>%
  ungroup() %>%
  filter(year %in% 2011:2020) %>%
  left_join(smoke_obs[,c("year", "fips","pop_2019", "death_rate_avg")])

month_bin <-  month_bin %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(month_bin))) %>%
  left_join(filter(coef_robust, model=="month bins")) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = nbins*(exp(coef)-1)*death_rate_avg/12*pop_2019) %>%
  group_by(year, bootid) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  group_by(bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% mutate(spec="monthly bins")


#--- annual bins in monthly conc
cutoff <- c(-0.01,0,0.25,0.5,0.75,1,2,3,4,5,10,20,150)
annual_month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  filter(age_group=="all_ages") %>%
  mutate(bins=cut(smokePM_mean, cutoff)) %>%
  ungroup() %>% group_by(fips, year, bins) %>%
  summarise(nbins=length(bins)) %>% ungroup()  %>%
  filter(year %in% 2011:2020) %>%
  left_join(smoke_obs[,c("year", "fips","pop_2019", "death_rate_avg")])

annual_month_bin <-  annual_month_bin %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(annual_month_bin))) %>%
  left_join(filter(coef_robust, model=="bins in monthly conc")) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = nbins*(exp(coef)-1)*death_rate_avg*pop_2019) %>%
  group_by(year, bootid) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  group_by(bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% mutate(spec="bins in monthly conc")

#--- annual bins in daily conc
cutoff <- c(-0.01,0.5,1,2,3,4,5,10,20,50, 2000)
annual_daily_bin <- readRDS(paste0(data_path, "/smokePM2pt5_predictions_daily_county_20060101-20201231.rds"))  %>%
  mutate(year=as.numeric(substr(date,1,4))) %>%
  mutate(fips=as.numeric(GEOID)) %>%
  mutate(bins=cut(smokePM_pred, cutoff)) %>% ungroup() %>%
  group_by(fips, year, bins) %>%
  summarise(nbins=length(bins)) %>% ungroup()  %>%
  filter(year %in% 2011:2020) %>%
  left_join(smoke_obs[,c("year", "fips","pop_2019", "death_rate_avg")]) %>%
  filter(!is.na(pop_2019))

annual_daily_bin <-  annual_daily_bin %>%
  slice(rep(1:n(), each = nboot)) %>%
  mutate(bootid=rep(1:nboot, nrow(annual_daily_bin))) %>%
  left_join(filter(coef_robust, model=="bins in daily conc")) %>%
  mutate(coef=replace(coef, is.na(coef), 0)) %>%
  mutate(death_smoke = nbins*(exp(coef)-1)*death_rate_avg*pop_2019) %>%
  group_by(year, bootid) %>% 
  summarise(death_smoke=sum(death_smoke)) %>% ungroup() %>%
  group_by(bootid) %>%
  summarise(death_smoke=mean(death_smoke)) %>% mutate(spec="bins in daily conc")

death_robust <- bind_rows(annual_daily_bin, annual_month_bin, month_bin, include_2020) %>%
  group_by(spec) %>%
  summarise(death_median=median(death_smoke),  
            death_mean=mean(death_smoke),
            death_p025=quantile(death_smoke, 0.025),
            death_p975=quantile(death_smoke, 0.975),
            death_p10=quantile(death_smoke, 0.1),
            death_p90=quantile(death_smoke, 0.9),
            death_p25=quantile(death_smoke, 0.25),
            death_p75=quantile(death_smoke, 0.75)) %>% ungroup() %>%
  bind_rows(main) %>%
  bind_rows(alternative_bins)

death_robust$spec <- factor(death_robust$spec, levels=c("main","alternative bins","include 2020","bins in daily conc", "bins in monthly conc","monthly bins"))
saveRDS(death_robust, paste0(health_path, "/summ_smoke_death_poisson_bin_robust.rds"))

death_robust <- readRDS(paste0(health_path, "/summ_smoke_death_poisson_bin_robust.rds")) %>% filter(spec!="monthly bins")

ggplot(death_robust, aes(x=spec, fill=spec, y=death_mean)) +
geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=death_p025, ymax=death_p975), stat = "identity",
                position = position_dodge2(.9, padding = .9), linewidth=0.9, width=0.01) +
  scale_fill_manual(values=c("grey40",met.brewer(name="Egypt", n=11, type="continuous")[c(1,5,7,9)]))+
  scale_colour_manual(values=c("grey40",met.brewer(name="Egypt", n=11, type="continuous")[c(1,5,7,9)]))+
  theme_classic() +
  theme(text = element_text(size=20)) +
  labs(x="", y = "Annual smoke deaths\n(2011-2020)", fill="") +
  theme(text = element_text(size=20), axis.text.x = element_blank()) 
ggsave("smoke_deaths_poisson_bin_robust.pdf", width=9, height=4.5)


###------------------------
## Archive: results for month models 
###-----------------------
cutoff <- c(-0.01,0,1,2,3,4,5,10,20,150)
month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

cutoff <- c(-0.01,0,1,2,5,10,150)
month_bin_coarse <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()
## monthly bins
fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          month_fe))

age_group <- unique(month_bin$age_group)
for (aa in age_group){
  df <- filter(month_bin, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:8,]) %>% 
    mutate(age_group=aa, panel="month", model="smoke conc bins", bins=levels(month_bin$smoke_bins)[2:9])
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          ma_fe))

age_group <- unique(month_bin$age_group)
for (aa in age_group){
  df <- filter(month_bin, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:8,]) %>% 
    mutate(age_group=aa, panel="month", model="smoke conc bins, Ma FE", bins=levels(month_bin$smoke_bins)[2:9])
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

## monthly bins (coarse)
fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          month_fe))

age_group <- unique(month_bin$age_group)
for (aa in age_group){
  df <- filter(month_bin_coarse, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:5,]) %>% 
    mutate(age_group=aa, panel="month", model="smoke conc bins (coarse)", bins=levels(month_bin_coarse$smoke_bins)[2:6])
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          ma_fe))

age_group <- unique(month_bin$age_group)
for (aa in age_group){
  df <- filter(month_bin_coarse, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:5,]) %>% 
    mutate(age_group=aa, panel="month", model="smoke conc bins, Ma FE (coarse)", bins=levels(month_bin_coarse$smoke_bins)[2:6])
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

coef_plot <- coef_combined %>% filter(panel=="month")
coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))
colnames(coef_plot)[1:2] <- c("beta","se")
ggplot(coef_plot, aes(x=age_group, colour=model)) +
  geom_point(aes(x=age_group, y=beta), position=position_dodge(width=0.2)) +
  geom_errorbar(aes(ymin=beta-1.96*se,
                    ymax=beta+1.96*se), position=position_dodge(width=0.2), width=0.01) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(x="", y="mortalities per 100k\nper 1ug smoke PM2.5") +
  geom_hline(aes(yintercept=0), linetype="dashed")
ggsave("county_month_mort.png", width=7, height=4)

#-------------------------------------------------------------------------
#----------------- Running log outcome models on annual outcomes ---------------
#-------------------------------------------------------------------------
mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))
#mort_smoke_month <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds"))

annual_fe <- "fips + state_name^year"

coef_combined <- NULL
age_group <- unique(mort_smoke_year$age_group)

### annual poissopn model
fmla <- as.formula(paste0("log(rate) ~ smokePM_mean + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

for (aa in age_group){
  df <- filter(mort_smoke_year, age_group == aa, year!=2020)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(t(model_tmp$coeftable[1,])) %>%
    mutate(age_group=aa, panel="year", model="log, linear")
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

#### poisson model on annual conc bins
cutoff <- c(-0.01,0.1,0.5,1,2,5,20)
annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

### annual bins
fmla <- as.formula(paste0("log(rate) ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

age_group <- unique(annual_bin$age_group)
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:5,]) %>% 
    mutate(age_group=aa, panel="year", model="log, smoke conc bin", bins=levels(annual_bin$smoke_bins)[2:6])
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}
write.xlsx(coef_combined, paste0(result_path, "/coef_annual_log.xlsx"))



#----------------- Running distributed lag model ----------------
month_lag <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(date=as.Date(paste(year, month, "01",sep="-"))) %>%
  group_by(death_type, fips, age_group) %>%
  arrange(date) %>%
  mutate(smoke_lag1=lag(smokePM_mean, n=1),
         smoke_lag2=lag(smokePM_mean, n=2),
         smoke_lag3=lag(smokePM_mean, n=3),
         smoke_lag4=lag(smokePM_mean, n=4),
         smoke_lag5=lag(smokePM_mean, n=5),
         smoke_lag6=lag(smokePM_mean, n=6),
         smoke_lag7=lag(smokePM_mean, n=7),
         smoke_lag8=lag(smokePM_mean, n=8),
         smoke_lag9=lag(smokePM_mean, n=9),
         smoke_lag10=lag(smokePM_mean, n=10),
         smoke_lag11=lag(smokePM_mean, n=11)) %>%
  ungroup()

filter(month_lag, fips==1001, age_group=="65_and_up") %>% 
  as.data.frame()

# fe="fips^month + fips^year"
# fe="state_name^month + state_name^year + fips"
fe="fips^month + state_name^year + fips + year^month"
fmla <- as.formula(paste0("rate ~ smokePM_mean + smoke_lag1 + smoke_lag2 +
                                  smoke_lag3 + smoke_lag4 + smoke_lag5 + smoke_lag6 +
                                  smoke_lag7 + smoke_lag8 + smoke_lag9 + smoke_lag10 + smoke_lag11 + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",fe))

coef_combined <- NULL
age_group <- unique(month_lag$age_group)
for (aa in age_group){
  df <- filter(month_lag, death_type == dd, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:12,]) %>% 
    mutate(death_type=dd, age_group=aa, panel="month", model="linear", lag=0:11)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

colnames(coef_combined)[1:2] <- c("beta", "se")

ggplot(coef_combined, aes(x=lag, y=beta)) +
  geom_point()+
  geom_errorbar(aes(ymin=beta-1.96*se,
                    ymax=beta+1.96*se),
                width=0.001) +
  facet_wrap(~age_group, ncol=1, scales="free") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic() +
  theme(text = element_text(size=14))
ggsave("month_lag_mortality_panel_MaFE.png", width=6, height=6)




