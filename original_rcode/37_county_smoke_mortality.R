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
# analysis of county smoke-mortality RELATIONSHIP
# Written by Minghao
# Last edited Dec 2023
#-------------------------------------------------------------------------------
mort_data <- readRDS(paste0(db_path, "/adaptation/data/mortality/processed/usa/age_standardized_rates/age_standardized_rates-county_month-all_sex-all_marital-all_cause.RDS")) %>%
  filter(year>2005, race_eth=="all") %>%
  rename(fips = fipsihme) %>%
  mutate(fips = as.numeric(fips))

mort_temp_data <-  readRDS(paste0(data_path, "/mortality_pwtemperature_panel_complete-fips.rds")) %>%
  filter(year>2005) %>%
  rename(fips = fipsihme) %>%
  mutate(fips = as.numeric(fips)) %>%
  distinct(year, month, fips, .keep_all = T)

met_list <- c("pw_temp_below_neg10", "pw_temp_neg10_neg5", "pw_temp_neg5_0", "pw_temp_0_5", "pw_temp_5_10", "pw_temp_10_15",
              "pw_temp_20_25", "pw_temp_25_30", "pw_temp_30_35", "pw_temp_above_35", "pwsum_prec")

##### combine mort data by age and temperature
mort_data_combined <- left_join(mort_data,
                                mort_temp_data[,c("year", "month", "fips", "state_name", met_list)],
                                by = c("year", "month", "fips")) %>%
  filter(!is.na(pw_temp_0_5))

county_smoke_month <- readRDS(paste0(data_path, "/smokePM2pt5_predictions_daily_county_20060101-20201231.rds")) %>%
           mutate(year=as.numeric(substr(date,1,4)),
                  month=as.numeric(substr(date,6,7))) %>%
  mutate(ndays=lubridate::days_in_month(date)) %>%
  group_by(GEOID, year, month, ndays) %>%
  summarise(smokePM_pred=sum(smokePM_pred)) %>% ungroup() %>%
  mutate(smokePM_mean=smokePM_pred/ndays) %>%
  mutate(fips=as.numeric(GEOID)) %>%
  select(-ndays, -smokePM_pred, -GEOID)

county_smoke_year <- readRDS(paste0(data_path, "/smokePM2pt5_predictions_daily_county_20060101-20201231.rds")) %>%
   mutate(year=as.numeric(substr(date,1,4))) %>%
  group_by(GEOID, year) %>%
  summarise(smokePM_pred=sum(smokePM_pred)) %>% ungroup() %>%
  mutate(n=if_else(leap_year(year),366,365),
         smokePM_mean=smokePM_pred/n) %>%
  mutate(fips=as.numeric(GEOID)) %>%
  select(-n, -smokePM_pred, -GEOID)

mort_smoke_month <- left_join(mort_data_combined, county_smoke_month, by=c("year", "month", "fips")) %>%
  mutate(smokePM_mean = replace(smokePM_mean, is.na(smokePM_mean), 0))

var_list <- c("rate",
              "pw_temp_below_neg10", "pw_temp_neg10_neg5", "pw_temp_neg5_0", "pw_temp_0_5", "pw_temp_5_10", "pw_temp_10_15",
              "pw_temp_20_25", "pw_temp_25_30", "pw_temp_30_35", "pw_temp_above_35", "pwsum_prec")

mort_year <- mort_data_combined %>% group_by(death_type, age_group, fips, year, pop, state_name) %>%
  summarise_at(var_list, .funs=sum, na.rm=T) %>%
  ungroup()

mort_smoke_year <- left_join(mort_year, county_smoke_year, by=c("year", "fips")) %>%
  mutate(smokePM_mean = replace(smokePM_mean, is.na(smokePM_mean), 0))

saveRDS(mort_smoke_year, paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))
saveRDS(mort_smoke_month, paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds"))

mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))
age_group <- unique(mort_smoke_year$age_group)

annual_fe <- "fips + state_name^year"
nboot <- 500

#------------------------------------------------
# linear model
#------------------------------------------------
fmla <- as.formula(paste0("rate ~ smokePM_mean + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
for (aa in age_group){
  df <- filter(mort_smoke_year, age_group == aa, year!=2020)
    for (bbb in 1:nboot){
        print(bbb)
        set.seed(bbb+100)

  df_tmp <- df[sample(1:nrow(df), replace = T), ]

  model_tmp <- feols(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
  coef_tmp <- data.frame(coef=model_tmp$coeftable[1,1]) %>%
    mutate(age_group=aa, panel="year", model="linear", bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
    }
}

saveRDS(coef_combined, paste0(result_path, "/coef_linear.rds"))

#------------ plot the coefficient 
# coef_combined <- read.xlsx(paste0(result_path, "/coef_county_linear_mort.xlsx")) %>%
#   filter(age_group!="under_5")
# 
# coef_plot <- coef_combined %>% filter(panel=="year")
# coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))
# colnames(coef_plot)[1:2] <- c("beta","se")
# ggplot(coef_plot, aes(x=age_group)) +
#   geom_point(aes(x=age_group, y=beta), position=position_dodge(width=0.2)) +
#   geom_errorbar(aes(ymin=beta-1.96*se,
#                     ymax=beta+1.96*se), position=position_dodge(width=0.2), width=0.01) +
#   theme_classic() +
#   theme(text = element_text(size=16)) +
#   labs(x="", y="mortalities per 100k\nper 1ug smoke PM2.5") +
#   geom_hline(aes(yintercept=0), linetype="dashed")
# ggsave("county_year_mort.pdf", width=5, height=3.5)

#-------------------------------------------------------------------------
#----------------- Running bin model to test non-linearity ---------------
#-------------------------------------------------------------------------
cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

cutoff <- c(-0.01,0.1, 0.25, 0.5, 0.75,1, 2.5, 5, 10, 20)
annual_bin_coarse <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

# annual_bin %>% distinct(fips, year, .keep_all = T) %>%
#   group_by(smoke_bins) %>%
#   summarise(n=length(year),
#             n_wo_2020=sum(as.numeric(year!=2020))) %>%
#   ungroup()

### annual bins (coarse)
fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin_coarse, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)

    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- feols(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:7,1]) %>% 
        mutate(age_group=aa, panel="year", model="linear bins (coarse)",
               bins=levels(annual_bin_coarse$smoke_bins)[2:8],
               bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

### annual bins (fine)
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)

    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- feols(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    coef_tmp <- data.frame(coef=model_tmp$coeftable[1:10,1]) %>% 
    mutate(age_group=aa, panel="year", model="linear bins", 
           bins=levels(annual_bin$smoke_bins)[2:11],
           bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

saveRDS(coef_combined, paste0(result_path, "/coef_linear_bins.rds"))


#------------- plotting coefs --------------
# coef_combined <- read.xlsx(paste0(result_path, "/coef_county_bins_mort.xlsx"))  %>%
#   filter(age_group!="under_5")
# 
# coef_plot <- coef_combined %>% filter(panel=="year", model=="smoke conc bin")
# coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))
# colnames(coef_plot)[1:2] <- c("beta","se")
# #coef_combined$bins <- factor(coef_combined$bins, levels=levels(month_bin$smoke_bins)[2:9])
# ggplot(coef_plot, aes(x=bins, y=beta)) +
#   geom_point()+
#   geom_errorbar(aes(ymin=beta-1.96*se,
#                     ymax=beta+1.96*se),
#                 width=0.001) +
#   facet_wrap(~age_group, ncol=1, scales="free") +
#   geom_hline(aes(yintercept=0), linetype="dashed") +
#   theme_classic() +
#   theme(text = element_text(size=14))
# ggsave("year_fine_bins_mortality_panel.pdf", width=5, height=6)
# 
# 
# coef_plot <- coef_combined %>% filter(panel=="month",
#                                       model %in% c("smoke conc bins", 
#                                                    "smoke conc bins, Ma FE"))
# coef_plot$age_group <- factor(coef_plot$age_group, levels= c("65_and_up", "under_65", "all_ages"))
# colnames(coef_plot)[1:2] <- c("beta","se")
# coef_plot$bins <- factor(coef_plot$bins, levels=unique(coef_plot$bins))
# ggplot(coef_plot, aes(x=bins, y=beta, colour=model)) +
#   geom_point(position=position_dodge(width=0.2))+
#   geom_errorbar(aes(ymin=beta-1.96*se,
#                     ymax=beta+1.96*se),
#                 width=0.001, 
#                 position=position_dodge(width=0.2)) +
#   facet_wrap(~age_group, ncol=1, scales="free") +
#   geom_hline(aes(yintercept=0), linetype="dashed") +
#   theme_classic() +
#   theme(text = element_text(size=14))
# ggsave("month_bins_mortality_panel.png", width=7.5, height=6)


#------------------------------------------------
#------- polynomial model in annual smoke
#------------------------------------------------
mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
mutate(smokePM_mean_sqr = smokePM_mean^2,
       smokePM_mean_cub = smokePM_mean^3)

### quadratic model
coef_combined <- NULL
fmla <- as.formula(paste0("rate ~ smokePM_mean + smokePM_mean_sqr + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

for (aa in age_group){
  df <- filter(mort_smoke_year, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)

  df_tmp <- df[sample(1:nrow(df), replace = T), ]
  
  model_tmp <- feols(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
  coef_tmp <- as.data.frame(t(model_tmp$coeftable[1:2,1])) %>%
    mutate(age_group=aa, panel="year", model="quadratic", bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}
saveRDS(coef_combined, paste0(result_path, "/coef_quadratic.rds"))


### cubic model
coef_combined <- NULL
fmla <- as.formula(paste0("rate ~ smokePM_mean + smokePM_mean_sqr + smokePM_mean_cub + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

for (aa in age_group){
  df <- filter(mort_smoke_year, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)

    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
    model_tmp <- feols(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    coef_tmp <- as.data.frame(t(model_tmp$coeftable[1:3,1])) %>%
      mutate(age_group=aa, panel="year", model="cubic", bootid=bbb)
    coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}
saveRDS(coef_combined, paste0(result_path, "/coef_cubic.rds"))


# coef_combined <- readRDS(paste0(result_path, "/coef_county_quadratic.rds"))
# #----- plot the response curve
# smoke_conc_vec <- seq(0,10, by=0.01)
# response <- NULL
# for (aa in age_group){
#   for (bbb in 1:nboot){
#     print(bbb)
#     coef_tmp <- filter(coef_combined, bootid==bbb, age_group==aa) 
#   response_tmp <- data.frame(smoke=smoke_conc_vec) %>%
#                    mutate(death_rate=as.numeric(coef_tmp[1])*smoke+
#                             as.numeric(coef_tmp[2])*smoke^2,
#                             #as.numeric(coef_tmp[3])*smoke^3,
#                           bottid=bbb, age_group=aa, model="quadratic")
#   response <- bind_rows(response_tmp, response)
#   }
#   }
# 
# response_summ <- response %>% group_by(age_group, smoke) %>%
# summarise(median=median(death_rate),
#           p025=quantile(death_rate, 0.025),
#           p975=quantile(death_rate, 0.975)) %>%
#   ungroup()
#   
# ggplot(response_summ %>% filter(smoke<9 ,age_group!="under_5")) +
#   geom_line(aes(x=smoke, y=median)) +
#   geom_ribbon(aes(x=smoke, ymin=p025, ymax=p975), alpha=0.5) +
#   facet_wrap(~age_group, ncol=1, scales="free") +
#   geom_hline(aes(yintercept=0), linetype="dashed") +
#   theme_classic() +
#   theme(text = element_text(size=18)) +
#   labs(x="", y="mortalities per 100k\nper 1ug smoke PM2.5")
# ggsave("year_quadratic.png", width=5, height=6)


#-------------------------------------------------------------------------
#----------------- Running poisson models on annual outcomes ---------------
#-------------------------------------------------------------------------
mort_smoke_year <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds"))

### annual poisson model
fmla <- as.formula(paste0("rate ~ smokePM_mean + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin_coarse, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
  model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)

  coef_tmp <- data.frame(coef=model_tmp$coeftable[1,1]) %>%
    mutate(age_group=aa, panel="year", model="poisson", bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}
saveRDS(coef_combined, paste0(result_path, "/coef_poisson.rds"))


#### poisson model on annual conc bins
cutoff <- c(-0.01, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6, 10, 20)
annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

cutoff <- c(-0.01,0.1, 0.25, 0.5, 0.75,1, 2.5, 5, 10, 20)
annual_bin_coarse <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

### annual bins
fmla <- as.formula(paste0("rate ~ i(smoke_bins)+ 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  for (bbb in 1:nboot){
    print(bbb)
    set.seed(bbb+100)
    
    df_tmp <- df[sample(1:nrow(df), replace = T), ]
    
  model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
    
  coef_tmp <- data.frame(coef=model_tmp$coeftable[1:10,1]) %>% 
    mutate(age_group=aa, panel="year", model="poisson bins", 
           bins=levels(annual_bin$smoke_bins)[2:11], bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
  }
}

### annual bins (coarse)
for (aa in age_group){
    df <- filter(annual_bin_coarse, age_group == aa, year!=2020)
 for (bbb in 1:nboot){
      print(bbb)
      set.seed(bbb+100)
      
      df_tmp <- df[sample(1:nrow(df), replace = T), ]
      
      model_tmp <- fepois(fml=fmla, weights = df_tmp$pop, cluster = df_tmp$fips, data = df_tmp)
      
  coef_tmp <- data.frame(coef=model_tmp$coeftable[1:7,1]) %>% 
    mutate(age_group=aa, panel="year", model="poisson bins (coarse)",
           bins=levels(annual_bin_coarse$smoke_bins)[2:8], bootid=bbb)
  coef_combined <- bind_rows(coef_combined, coef_tmp)
 }
}

saveRDS(coef_combined, paste0(result_path, "/coef_poisson_bins.rds"))



#-------------------------------------------------------------------------
#----------------- Running bin model in monthly conc to test non-linearity ---------------
#-------------------------------------------------------------------------
cutoff <- c(-0.01,0,1,2,3,4,5,10,20,150)
month_bin <- readRDS(paste0(data_path, "/county_mort_smoke_month_2006_2020_age.rds")) %>%
  mutate(smoke_bins=cut(smokePM_mean, cutoff)) %>%
  ungroup()

annual_month_bin <- month_bin %>% group_by(fips, year, smoke_bins, age_group) %>%
  summarise(n=length(smoke_bins)) %>% ungroup() %>%
  pivot_wider(id_cols = fips:age_group, names_from = smoke_bins, values_from = n)
colnames(annual_month_bin)[4:12] <- c("smoke_0","smoke_0_1", "smoke_1_2",
                                      "smoke_4_5","smoke_5_10", "smoke_3_4",
                                      "smoke_2_3", "smoke_10_20", "smoke_20_plus")

annual_bin <- readRDS(paste0(data_path, "/county_mort_smoke_year_2006_2020_age.rds")) %>%
  left_join(annual_month_bin) %>%
  mutate(across(c("smoke_0","smoke_0_1", "smoke_1_2",
                  "smoke_4_5","smoke_5_10", "smoke_3_4",
                  "smoke_2_3", "smoke_10_20", "smoke_20_plus"), ~replace_na(.x, 0)))

fmla <- as.formula(paste0("rate ~ smoke_0_1 + smoke_1_2 +  smoke_2_3 + smoke_3_4 +
                                  smoke_4_5 +smoke_5_10 + smoke_10_20 + smoke_20_plus + 
                          pw_temp_below_neg10 + pw_temp_neg10_neg5 + pw_temp_neg5_0 + pw_temp_0_5 +pw_temp_5_10 + pw_temp_10_15 +
                          pw_temp_20_25 + pw_temp_25_30 + pw_temp_30_35 + pw_temp_above_35 + pwsum_prec + pwsum_prec^2 | ",
                          annual_fe))

coef_combined <- NULL
age_group <- unique(annual_bin$age_group)
for (aa in age_group){
  df <- filter(annual_bin, age_group == aa, year!=2020)
  print(aa)
  model_tmp <- feols(fml=fmla, weights = df$pop, cluster = df$fips, data = df)
  coef_tmp <- as.data.frame(model_tmp$coeftable[1:8,]) %>% 
    mutate(age_group=aa, panel="year", model="bin of smoke months")
  coef_tmp$bins <- row.names(coef_tmp)
  
  coef_combined <- bind_rows(coef_combined, coef_tmp)
}

write.xlsx(coef_combined, paste0(result_path, "/coef_county_annual_bins_smokemonth.xlsx"))


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








Poisson

Log(rate)

Fepois



