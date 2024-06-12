library(terra);library(ncdf4)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(plotly)
library(sp);library(sf)
library(dplyr)   
library(fixest)
library(tidyverse)
library(MetBrewer)
library(stringr)
library(splines2)
library(exactextractr)
library(zoo)

rm(list=ls())
gc()

code_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/rcode"
data_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/data"
figure_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/figures"
result_path <- "/Users/mhqiu/Dropbox_Personal/Dropbox/Postdoc/Research/smoke_climate/results/smoke_model"

db_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects/smoke-climate"
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"

setwd(figure_path)
states_map <- map_data("state")

pop_grid <- readRDS(paste0(paste0(data_path, "/health/us_population_grid_2020_2022_2050.rds")))
smokegrid_climate <- readRDS(paste0(data_path, "/crosswalk_10km_grid_US_9region.rds"))
pop_df <- as.data.frame(pop_grid) %>% select(grid_id_10km, pop_2020, pop_2022, pop_2050) %>%
  left_join(smokegrid_climate[,c("grid_id_10km", "NAME", "region")])

# color_scenario <- c("grey65",
#                     rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])
color_scenario <- c("grey65", "#1F2E52", "#E78928","#CB172F") # "#38A2C6", "#870A21"

#-------------------------------------------------------------------------------
# Plot the figure 3
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------

#---------------------------------------------
#--- Fig 3A: plot the smoke maps
#---------------------------------------------
# read the smoke grid shape file
grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) %>%
  rename(grid_id_10km=ID)

smoke_2050_df <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
  filter(year %in% c(2055)) %>% mutate(year=as.character(year))

smoke_2011_2020 <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>% 
  filter(year %in% 2011:2020) %>%
  group_by(grid_id_10km) %>%
  summarise_at(c("smokePM_mean"),.funs = mean, na.rm=T) %>%
  mutate(scenario="Observation", year="2011-2020") %>%
  rename(smoke_grid = smokePM_mean) %>% filter(grid_id_10km %in% unique(smoke_2050_df$grid_id_10km))

smoke_combined <- bind_rows(smoke_2050_df, smoke_2011_2020) %>%
  group_by(year, scenario, grid_id_10km) %>%
  summarise(smoke_grid = median(smoke_grid)) %>% ungroup()

#### Plotting ######
col_disc <- c("grey90","#c3f4f6","#fff179","#f5b642","orangered3","#551f00")

smoke_tmp <- right_join(grid_10km, smoke_combined)

label_smoke <- c(0,0.25,1,2.5,5,10, max(smoke_tmp$smoke_grid, na.rm=T))
smoke_tmp$smoke_cut <- cut(smoke_tmp$smoke_grid, label_smoke, include.lowest = T)
smoke_tmp$smoke_cut <- factor(smoke_tmp$smoke_cut, levels= rev(levels(smoke_tmp$smoke_cut)))

ggplot(smoke_tmp, aes(colour=smoke_cut, fill=smoke_cut)) +
  geom_sf() + 
  geom_polygon(data=states_map,aes(long, lat, group = group),
               fill=NA,color ="black",size=0.3) +
    scale_color_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  labs(colour="smoke (ug m-3)", fill="smoke (ug m-3)") +
  facet_wrap(~year+scenario)
ggsave("median_GCMs/proj_smokePM_map_gridmetobs.png", width=12, height=8)

##How unprecedented is condition of >1ug
data_2050 <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
  filter(year %in% c(2055)) %>% mutate(year=as.character(year)) %>%
  group_by(year, scenario, grid_id_10km) %>%
  summarise(smoke_grid = median(smoke_grid)) %>% ungroup()

grid_2050 <- filter(data_2050, scenario=="ssp370", smoke_grid>1)$grid_id_10km %>% unique()

data_2011_2020 <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds")) %>%
  filter(grid_id_10km %in%grid_2050, smokePM_mean>1, year!=2020)
grid_2011_2020  <- unique(data_2011_2020$grid_id_10km)

#---------------------------------------------
#--- Fig 3B: plot the smoke at the state-level
#---------------------------------------------
smoke_2050 <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
  filter(year %in% c(2055)) %>% mutate(year=as.character(year)) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020)) %>%
  rename(state=NAME)

smoke_2050_state <- smoke_2050 %>% group_by(year, scenario, gcm ,state) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020),  ## population-weighted does not matter since 2050 pop is just scaling up
            pop_total = weighted.mean(totalPM, w=pop_2020)) %>%
  ungroup() %>%
  mutate(smoke_percent=pop_smoke/pop_total) %>%
  group_by(year, scenario, state) %>%
  summarise(smoke_percent=median(smoke_percent),
            pop_smoke=median(pop_smoke),
            pop_total=median(pop_total)) %>% ungroup()

smoke_2011_2020 <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>% 
  filter(year %in% 2011:2020) %>%
  rename(smoke_grid = smokePM_mean) %>% 
  filter(grid_id_10km %in% unique(smoke_2050$grid_id_10km)) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020)) %>%
  rename(state=NAME)

smoke_2020_state <- smoke_2011_2020 %>% group_by(state) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020),  ## population-weighted does not matter since 2050 pop is just scaling up
            pop_total = weighted.mean(totalPM, w=pop_2020)) %>%
  ungroup() %>%
  mutate(smoke_percent=pop_smoke/pop_total)   %>%
  mutate(scenario="Observation", year="2011-2020")

smoke_pop_state <- bind_rows(smoke_2020_state, smoke_2050_state)


state_list_10 <- smoke_pop_state %>% filter(scenario=="ssp370") %>%
  slice_max(smoke_percent, n=10) %>% select(state) %>% unlist()

ggplot(smoke_pop_state %>% filter(state %in% state_list_10)) +
  geom_bar(aes(x=reorder(state, smoke_percent), y=smoke_percent, fill=scenario), 
           stat="identity", position=position_dodge()) +
  coord_flip() +
  theme_classic() + 
  theme(text = element_text(size=21),
        axis.text = element_text(size=21)) +
  labs(x="", y="Smoke contributions to total pop-weighted PM") +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values=color_scenario)
  ggsave("median_GCMs/state_top10_smoke_percent.pdf", width=8, height=6)

#### output information for all the states
state_list <- smoke_pop_state %>% filter(scenario=="ssp370") %>%
  arrange(desc(smoke_percent)) %>% select(state) %>% unlist()
smoke_pop_state$state <- factor(smoke_pop_state$state, levels = state_list)
smoke_pop_state <- smoke_pop_state %>% arrange(state,scenario) %>%
  mutate(scenario=replace(scenario, scenario=="Observation", "2011-2020"),
         scenario=replace(scenario, scenario=="ssp126", "SSP1-2.6"),
         scenario=replace(scenario, scenario=="ssp245", "SSP2-4.5"),
         scenario=replace(scenario, scenario=="ssp370", "SSP3-7.0"))
write.xlsx(smoke_pop_state, paste0(result_path, "/smoke_pecent_state_median.xlsx"))

#---------------------------------------------
#--- Fig 3C: plot the population-weighted smoke
#---------------------------------------------
  smoke_future <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
    left_join(pop_df) %>%
    filter(!is.na(smoke_grid), !is.na(pop_2020)) 
  
  smoke_future_us <- smoke_future %>% group_by(year, scenario, gcm) %>%
    summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020)) %>%
    ungroup() %>%
    group_by(year, scenario) %>% 
    summarise(pop_smoke=median(pop_smoke)) %>% ungroup() %>% ## mean or median across GCMs
    mutate(year=year-5)
  
  smoke_2011_2020 <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>% 
    filter(year %in% 2011:2020) %>%
    rename(smoke_grid = smokePM_mean) %>% 
    filter(grid_id_10km %in% unique(smoke_future$grid_id_10km)) %>%
    left_join(pop_df) %>%
    filter(!is.na(smoke_grid), !is.na(pop_2020))
  
  smoke_2020_us <- smoke_2011_2020 %>%
    summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020)) %>%
    ungroup() %>%
    mutate(scenario="Observation", year=2020)

### just for visualization 
smoke_pop_us <- bind_rows(smoke_future_us %>% filter(year!=2020), 
                          smoke_2020_us,
                          smoke_2020_us %>% mutate(scenario="ssp126"),
                          smoke_2020_us %>% mutate(scenario="ssp245"),
                          smoke_2020_us %>% mutate(scenario="ssp370"))

ggplot() +
  geom_point(data = smoke_pop_us, 
             aes(x=year, y=pop_smoke, colour=scenario), size=3, position = position_dodge(width=0.1)) +
  geom_line(data = smoke_pop_us %>% filter(scenario %in% c("ssp126", "ssp245", "ssp370")), 
            aes(x=year, y=pop_smoke, colour=scenario), size=1.4) +
   geom_point(data = smoke_pop_us %>% filter(scenario == "ssp126", year==2020), 
             aes(x=year, y=pop_smoke), colour="grey65", size=3) +
   scale_colour_manual(values=color_scenario) +
  scale_x_continuous(breaks=c(2020,2030,2040,2050),
                     labels=c("2011-2020",2030,2040,2050)) +
  theme_classic() +
  theme(text = element_text(size=18),
        axis.text = element_text(size=16)) +
  scale_y_continuous(breaks=c(0.5, 1, 1.5), limits = c(0.2,1.5)) +
  labs(x= "", y="Pop-weighted smoke PM2.5")
ggsave("median_GCMs/pop_smoke_us_timeseries.pdf", width=6.5, height=3.5)

#---------------------------------------------
#--- Fig 3d: Uncertainty over GCMs
#---------------------------------------------
smoke_future <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020)) %>%
  filter(year==2055)

smoke_future_us <- smoke_future %>% group_by(year, scenario, gcm) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020)) %>%
  ungroup() 

smoke_pop_gcm_quantile <- smoke_future_us %>% group_by(year, scenario) %>%
  summarise(smoke_50=median(pop_smoke),
            smoke_mean=mean(pop_smoke),
            smoke_10=quantile(pop_smoke, 0.1),
            smoke_25=quantile(pop_smoke, 0.25),
            smoke_75=quantile(pop_smoke, 0.75),
            smoke_90=quantile(pop_smoke, 0.9)) %>%
  ungroup()

ggplot(smoke_pop_gcm_quantile, aes(x=scenario)) +
  geom_errorbar(aes(ymin=smoke_25, ymax=smoke_75,colour=scenario), size=5, width=0.0001) +
  geom_errorbar(aes(ymin=smoke_10, ymax=smoke_90,colour=scenario), width=0.0001) +
 # geom_errorbar(aes(ymin=smoke_50, ymax=smoke_50), colour="black",width=0.55) +
  scale_colour_manual(values=color_scenario[2:4]) +
  geom_point(aes(y=smoke_50), colour="yellow",size=2.5) +
  #scale_fill_manual(values=color_scenario) +
  theme_classic() +
  theme(text = element_text(size=18),
        axis.text = element_text(size=18),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(x="") +
  scale_y_continuous(breaks=seq(0.5,2.5,by=0.5), limits=c(0.4,3.1))
ggsave("median_GCMs/2050_10yrs_gcm_uncertainty_25_10.pdf", width=3.1,height=4)


#---------------------------------------------
# Fig 3E: plot the variability over individual years
#---------------------------------------------
smoke_2046_2055<- readRDS(paste0(result_path, "/smoke_proj/smoke_2046_2055_debias_gcm.rds")) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020))

smoke_pop_2046_2055 <- smoke_2046_2055 %>% group_by(gcm, year, scenario) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2050)) %>% ungroup() %>%
  group_by(year, scenario) %>%
  summarise(pop_smoke_median=median(pop_smoke),
            pop_smoke_mean=mean(pop_smoke)) %>% ungroup()

smoke_pop_2046_2055 %>% group_by(scenario) %>%
  summarise(mean=mean(pop_smoke_mean),
            median=mean(pop_smoke_median)) %>% ungroup()

ggplot() +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp126"),
                 aes(y=pop_smoke_mean, x=0, xmin=-1, xmax=1, colour=scenario), 
                 size=1, height=0.001) +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp245"),
                 aes(y=pop_smoke_mean, x=4, xmin=3, xmax=5,colour=scenario), size=1, height=0.001) +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp370"),
                 aes(y=pop_smoke_mean, x=8, xmin=7, xmax=9,colour=scenario), size=1, height=0.001) +
  scale_colour_manual(values=color_scenario[2:4]) +
  labs(x="",y="") + 
  theme_classic() +
  theme(text = element_text(size=18),
        axis.text = element_text(size=18),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_y_continuous(breaks=seq(0.5,2.5,by=0.5), limits=c(0.4,3.1))
ggsave("2046_2055_years_uncertainty.pdf", width=3.6,height=4)

###calculate the relative differences across GCMs
smoke_2046_2055<- readRDS(paste0(result_path, "/smoke_proj/smoke_2046_2055_debias_gcm.rds")) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020))

smoke_pop_2046_2055 <- smoke_2046_2055 %>% group_by(gcm, year, scenario) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2050)) %>% ungroup() %>%
  group_by(gcm, scenario) %>%
  arrange(pop_smoke) %>%
  mutate(index= row_number()) %>% ungroup() %>% ## for each GCM and year creat index sorted by pop_smoke
  group_by(index, scenario) %>%
  summarise(pop_smoke_median=median(pop_smoke),
            pop_smoke_mean=mean(pop_smoke)) %>% ungroup()

ggplot() +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp126"),
                 aes(y=pop_smoke_median, x=0, xmin=-1, xmax=1, colour=scenario), 
                 size=1, height=0.001) +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp245"),
                 aes(y=pop_smoke_median, x=4, xmin=3, xmax=5,colour=scenario), size=1, height=0.001) +
  geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp370"),
                 aes(y=pop_smoke_median, x=8, xmin=7, xmax=9,colour=scenario), size=1, height=0.001) +
  scale_colour_manual(values=color_scenario[2:4]) +
  labs(x="",y="") + 
  theme_classic() +
  theme(text = element_text(size=18),
        axis.text = element_text(size=18),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  scale_y_continuous(breaks=seq(0.5,4.5,by=1), limits=c(0.4,5.1))
ggsave("median_GCMs/2046_2055_years_uncertainty_byindex.pdf", width=3.6,height=4)


### calculate the ratios between bad/good year
smoke_pop_2046_2055 <- smoke_2046_2055 %>% group_by(gcm, year, scenario) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2050)) %>% ungroup() %>%
  group_by(gcm, scenario) %>%
  mutate(ratio=max(pop_smoke)/min(pop_smoke)) %>% ungroup() 

ratio_gcm_quantile <- smoke_pop_2046_2055  %>% group_by(year, scenario) %>%
  summarise(ratio_50=median(ratio),
            ratio_mean=mean(ratio),
            ratio_10=quantile(ratio, 0.1),
            ratio_25=quantile(ratio, 0.25),
            ratio_75=quantile(ratio, 0.75),
            ratio_90=quantile(ratio, 0.9)) %>%
  ungroup()

ggplot(ratio_gcm_quantile, aes(x=scenario)) +
  # geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp126"),
  #                aes(y=ratio, x=0, xmin=-1, xmax=1, colour=scenario), 
  #                size=1, height=0.001) +
  # geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp245"),
  #                aes(y=ratio, x=4, xmin=3, xmax=5,colour=scenario), size=1, height=0.001) +
  # geom_errorbarh(data=smoke_pop_2046_2055 %>% filter(scenario=="ssp370"),
  #                aes(y=ratio, x=8, xmin=7, xmax=9,colour=scenario), size=1, height=0.001) +
  geom_errorbar(aes(ymin=ratio_25, ymax=ratio_75,colour=scenario), size=5, width=0.0001) +
  geom_errorbar(aes(ymin=ratio_10, ymax=ratio_90,colour=scenario), width=0.0001) +
  # geom_errorbar(aes(ymin=ratio_50, ymax=ratio_50), colour="black",width=0.55) +
  scale_colour_manual(values=color_scenario[2:4]) +
  geom_point(aes(y=ratio_50), colour="yellow",size=2.5) +
  scale_colour_manual(values=color_scenario[2:4]) +
  labs(x="",y="") + 
  theme_classic() +
  theme(text = element_text(size=18),
        axis.text = element_text(size=18),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 
ggsave("2046_2055_years_uncertainty_ratio.pdf", width=3.6,height=4)

#### calculate the smoke contribution in total pop-weighted PM in today and future
smoke_2050 <- readRDS(paste0(result_path, "/smoke_proj/smoke_10yrs_debias_gcm.rds")) %>% 
  filter(year %in% c(2055)) %>% mutate(year=as.character(year)) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020)) 

smoke_2050_us <- smoke_2050 %>% group_by(year, scenario, gcm) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020),  ## population-weighted does not matter since 2050 pop is just scaling up
            pop_total = weighted.mean(totalPM, w=pop_2020)) %>%
  ungroup() %>%
  mutate(smoke_percent=pop_smoke/pop_total) %>%
  group_by(year, scenario) %>%
  summarise(smoke_percent=median(smoke_percent)) %>% ungroup()

smoke_2011_2020 <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>% 
  filter(year %in% 2011:2020) %>%
  rename(smoke_grid = smokePM_mean) %>% 
  filter(grid_id_10km %in% unique(smoke_2050$grid_id_10km)) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020)) %>%
  rename(state=NAME)

smoke_2020_us <- smoke_2011_2020 %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020),  ## population-weighted does not matter since 2050 pop is just scaling up
            pop_total = weighted.mean(totalPM, w=pop_2020)) %>%
  ungroup() %>%
  mutate(smoke_percent=pop_smoke/pop_total)   %>%
  mutate(scenario="Observation", year="2011-2020")

smoke_pop_us <- bind_rows(smoke_2020_us, smoke_2050_us)

#### calculate pop-weight effects by year 
smoke_obs <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2006_2020.rds"))  %>% 
  rename(smoke_grid = smokePM_mean) %>% 
  filter(grid_id_10km %in% unique(smoke_2050$grid_id_10km)) %>%
  left_join(pop_df) %>%
  filter(!is.na(smoke_grid), !is.na(pop_2020))

smoke_2020_us <- smoke_obs %>%
  group_by(year) %>%
  summarise(pop_smoke = weighted.mean(smoke_grid, w=pop_2020)) %>%
  ungroup()


####### Smoke contribution at grid cell level ###########

col_disc <- c("grey85","#c3f4f6","#fff179","#f5b642","orangered3","#551f00")
smoke_tmp <- right_join(grid_10km, smoke_2050_df)

label_smoke <- c(0,0.05,0.1,0.2,0.5,0.75, 1)
smoke_tmp$smoke_cut <- cut(smoke_tmp$smoke_percent, label_smoke, include.lowest = T)
smoke_tmp$smoke_cut <- factor(smoke_tmp$smoke_cut, levels= rev(levels(smoke_tmp$smoke_cut)))

col_map <- c(rep("white",2), rev(met.brewer(name="Demuth", n=10, type="continuous")[1:6]))
ggplot() +
  geom_sf(data=smoke_tmp %>% filter(year %in% c("2016-2020", "2055")),
          aes(colour=smoke_cut, fill=smoke_cut), linewidth=0.0001) +
    geom_polygon(data=states_map,aes(long, lat, group = group),
               fill=NA,color ="black",size=0.3) +
  scale_colour_manual(values=rev(col_disc)) +
  scale_fill_manual(values=rev(col_disc)) +
   theme_classic() + theme(text = element_text(size=14),
                          axis.ticks = element_blank(),
                          axis.line = element_blank(),
                          axis.text = element_blank()) +
  facet_wrap(~year+scenario)
ggsave("smoke_percent_proj_gridmetobs.png", width=12, height=8)



#####--------- ARCHIVE
### read observed smoke between 2016 and 2020
# total_smoke_pm <- readRDS(paste0(data_path, "/totalPM/totalPM_smoke_combined_2016_2020.rds")) 
# smoke_2016_2020 <- total_smoke_pm %>% 
#   group_by(grid_id_10km) %>%
#   summarise_at(c("totalPM","smokePM_mean"),.funs = mean, na.rm=T) %>%
#   ungroup() %>%
#   mutate(nonsmoke = totalPM - smokePM_mean)
# smoke_2016_2020[smoke_2016_2020$nonsmoke<0, "nonsmoke"] <- 0
# 
# ### read the projected smoke in 2055
# smoke_2050 <- list.files(path=paste0(result_path, "/smoke_proj/mean_10yrs/"),
#                          pattern="proj_smoke_gridmetobs_2050_10yrs*",
#                          full.names = T) %>%
#   map(readRDS) %>%
#   bind_rows() %>%
#   group_by(grid_id_10km, year, scenario, gcm) %>%
#   summarise(pred_smoke=mean(pred_smoke)) %>% ungroup() %>%
#   group_by(grid_id_10km, year, scenario) %>%
#   summarise(pred_smoke=median(pred_smoke)) %>% ungroup()
# 
# 
# ### read the predicted smoke using the same reg model  
# pred_smoke_gfed <- readRDS(paste0(result_path, "/smoke_proj/gfedDM_final_smoke_pred_grid_wind_9region_90cone_gridmet.rds")) %>%
#   filter(year %in% 2016:2020) %>%
#   filter(model %in% c("wind_all", "obs")) %>%
#   group_by(grid_id_10km, model) %>%
#   summarise(pred_smoke=mean(pred_smoke, na.rm=T)) %>%
#   ungroup() %>%
#   pivot_wider(values_from=pred_smoke, names_from = model, names_prefix = "smoke_") 
# 
# ##de-bias based on monthly-mean in 2016-2020
# smoke_2050_debias <- left_join(smoke_2050, pred_smoke_gfed) %>%
#   mutate(pred_smoke_debias=pred_smoke - smoke_wind_all + smoke_obs) %>%
#   rename(smoke_grid = pred_smoke_debias ) %>%
#   mutate(year=as.character(year))
# 
# grid_list <- intersect(unique(smoke_2050_debias$grid_id_10km), unique(smoke_2016_2020$grid_id_10km))
# 
# smoke_2016_2020 <- filter(smoke_2016_2020, grid_id_10km%in%grid_list) %>%
#   mutate(scenario="Observation",
#          year="2016-2020") %>%
#   rename(smoke_grid = smokePM_mean)
# 
# smoke_2050_df <- bind_rows(smoke_2050_debias, smoke_2016_2020) %>% filter(!is.na(smoke_grid)) %>%
#   filter(grid_id_10km %in% grid_list) %>% select(grid_id_10km,  year, scenario, smoke_grid)
# smoke_2050_df[smoke_2050_df$smoke_grid<0, "smoke_grid"] <- 0
# 
# smoke_2050_df <- left_join(smoke_2050_df, smoke_2016_2020[,c("grid_id_10km", "nonsmoke")]) %>%
#   mutate(totalPM=smoke_grid+nonsmoke,
#          smoke_percent=smoke_grid/totalPM)
# saveRDS(smoke_2050_df, paste0(result_path, "/smoke_proj/smoke_2050_10yrs_debias_percent.rds"))

