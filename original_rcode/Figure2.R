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

setwd(figure_path)
states_map <- map_data("state")

color_map <- c("grey65",
               rev(met.brewer(name="Homer1", n=8, type="discrete"))[c(2,6,8)])

#-------------------------------------------------------------------------------
# Plot the figure 2
# Written by Minghao
# Last edited Jan 2024
#-------------------------------------------------------------------------------


#---------------------------------------------
#--- Fig 2A: plot the smoke model coefficients
#---------------------------------------------
coef_region <- readRDS(paste0(result_path,"/smoke_proj/gfedDM_final_coef_wind_9region_90cone_gridmet.rds"))  %>% 
  filter(model=="wind_all", 
         var!="year") %>%
  mutate(dist=gsub("^([^.]+)km.*", "\\1", var)) %>% 
  mutate(dist=gsub(".*dist([^.]+)$", "\\1", dist)) %>% 
  mutate(dist=gsub(".*above_([^.]+)$", "\\1", dist)) %>%
  mutate(wind=gsub(".*km_([^.]+)$", "\\1", var)) %>% 
  mutate(wind=replace(wind, !wind%in%c("up","other","down"), "nowind"))

bin_dist<- c("within50", "50_100", "100_200", "200_350","350_500","500_750",  
             "750_1000","1000_1500","1500_2000", "2000")
label_dist <- c("<50km", "50-100", "100-200", "200-350","350-500","500-750",  
                "750-1000","1000-1500","1500-2000", ">2000km")
coef_region$dist <- factor(coef_region$dist, levels = bin_dist, labels = label_dist)
coef_region <- coef_region %>% mutate(se=se*1e8, beta=beta*1e8)  ###for GFED multiply by 1e8 

region_list <- c("Northwest","West","Northern Rockies","Southwest","South","Upper midwest","Ohio valley","Southeast","Northeast" )
coef_region$region <- factor(coef_region$region, levels=region_list)

#col_map <-  met.brewer(name="Manet", n=12, type="continuous")[c(8,4,2)]
col_map <-  met.brewer(name="VanGogh1", n=7, type="continuous")[c(7,5,1)]
ggplot(coef_region %>% filter(region=="Northern Rockies"), 
       aes(x=dist, y=beta, colour=wind)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), 
                position=position_dodge(width=0.5), width=0, size=1) +
   theme_classic() + 
  theme(text=element_text(size=20), axis.text.x = element_text(angle=90)) +
  labs(x="",y="ug smoke per\n100,000 ton DM") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  scale_color_manual(values = col_map) +
  facet_wrap(~region)
ggsave("northrocky_smoke_coef.pdf",width=8.5, height = 5)

ggplot(coef_region %>% filter(
  region=="Northern Rockies",
  !dist %in% c("<50km", "50-100", "100-200", "200-350")), 
       aes(x=dist, y=beta, colour=wind)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), 
                position=position_dodge(width=0.5), width=0, size=1) +
  theme_classic() + 
  theme(text=element_text(size=20), axis.text.x = element_text(angle=90)) +
  labs(x="",y="ug smoke per\n100,000 ton DM") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  scale_color_manual(values = col_map) +
  facet_wrap(~region)
ggsave("northrocky_smoke_coef_longdist.pdf",width=7.5, height = 5)


#---------------------------------------------------------------------
#--- Fig 2B: show the regional heterogeneity of the smoke model coefs
#---------------------------------------------------------------------
col_9region <-  met.brewer(name="Juarez", n=9, type="continuous")

ggplot(coef_region %>% filter(
  dist %in% c("<50km", "50-100"),
  wind == "up"), 
  aes(x=wind, y=beta, colour=region)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), 
                position=position_dodge(width=0.5), width=0, size=1) +
  theme_classic() + 
  theme(text=element_text(size=22), axis.text = element_text(size=22)) +
  labs(x="",y="ug smoke per\n100,000 ton DM") + 
  geom_hline(aes(yintercept=0), linetype="dashed") +
  scale_color_manual(values = col_9region) +
  facet_wrap(~dist, scales = "free",nrow=1) +
  guides(colour=guide_legend(ncol=2)) +
  ylim(-0.2,2.3)
ggsave("region_het_50_100.pdf", width=12.5 ,height = 4)


### plot the map of the 9 regions
smokegrid_climate <- readRDS(paste0(data_path, "/crosswalk_10km_grid_US_9region.rds"))
db_proj_path <- "/Users/mhqiu/BurkeLab Dropbox/Projects"
grid_10km <- st_read(paste0(db_proj_path, "/smoke_PM_prediction/data/1_grids/grid_10km_wgs84/grid_10km_wgs84.shp")) %>%
  rename(grid_id_10km=ID)
grid_10km_region <- inner_join(grid_10km, smokegrid_climate[c("grid_id_10km", "region")])
grid_10km_region$region <- factor(grid_10km_region$region, levels=region_list)

ggplot(grid_10km_region) +
  geom_sf(aes(fill=region, colour=region), size=0.001) +
  theme_void() +
  scale_colour_manual(values=col_9region) +
  scale_fill_manual(values=col_9region) 
ggsave("climate_9_regions.png", width=6, height=3)


#Figure Sx: plot the R2
smoke_r2 <-  readRDS(paste0(result_path,"/smoke_proj/gfedDM_final_r2_wind_9region_90cone_gridmet.rds")) %>%
  filter(model == "wind_all")

ggplot(smoke_r2) +
  geom_bar( aes(x=reorder(region, withinr2), y=withinr2, fill=region), stat = "identity") +
  geom_point(aes(x=reorder(region, adjr2), y=adjr2), colour="black", size=2) +
  theme_classic() + 
  theme(text=element_text(size=23), axis.text.x = element_text(angle=90)) +
  labs(x="", y="R-Squared") + 
  scale_y_continuous(limits=c(0,0.85)) +
  scale_fill_manual(values=col_9region)
ggsave("gfedPM25_region9_wind_R2.png",width=11, height = 6)


# ################################  reg results by 9 regions   ################################
# r2_9region <- readRDS(paste0(result_path,"/gfedPM25_r2_wind_9region_90cone.rds")) 
# r2_9region$model <- factor(r2_9region$model,  levels = c("nomet","nowind","wind_main","wind_coarse","wind_all","wind_all01"),
#                            labels = c("No met, no wind","No wind","By wind","By wind (coarse)","By wind (all)","By wind (only up+down)"))
# 
# ggplot(filter(r2_9region, model%in%c("By wind (all)"), region!="US")) +
#   geom_bar( aes(x=reorder(region, withinr2), y=withinr2, fill=region), stat = "identity") +
#   geom_point(aes(x=reorder(region, adjr2), y=adjr2), colour="black", size=2) +
#   theme_bw() + theme(text=element_text(size=20), axis.text.x = element_text(angle=90)) +
#   labs(x="",y="Within R2") + geom_hline(aes(yintercept=0), linetype="dashed") + scale_y_continuous(limits=c(0,0.85))
# ggsave("gfedPM25_region9_wind_R2.png",width=11, height = 5)
# 
# coef_9region <- readRDS(paste0(result_path,"/gfedPM25_coef_wind_9region_90cone.rds")) 
# coef_9region$model <- factor(coef_9region$model,  levels = c("nomet","nowind","wind_main","wind_coarse","wind_all","wind_all01"),
#                              labels = c("No met, no wind","No wind","By wind","By wind (coarse)","By wind (all)","By wind (only up+down)"))
# 
# coef_9region <- coef_9region %>% ##bind_rows(coef_9region, coef_us)
#   filter(var!="year") %>%
#   mutate(dist=gsub("^([^.]+)km.*", "\\1", var)) %>% 
#   mutate(dist=gsub(".*dist([^.]+)$", "\\1", dist)) %>% 
#   mutate(dist=gsub(".*above_([^.]+)$", "\\1", dist)) %>%
#   #mutate(dist=str_sub(dist, 1, 10)) %>%
#   mutate(wind=gsub(".*km_([^.]+)$", "\\1", var)) %>% 
#   mutate(wind=replace(wind, !wind%in%c("up","other","down"), "nowind"))
# 
# bin_dist<- c("within50", "50_100", "100_200", "200_350","350_500","500_750",  
#              "750_1000","1000_1500","1500_2000", "2000",  "50_200", "200_500", "500_1000")
# label_dist <- c("<50km", "50-100", "100-200", "200-350","350-500","500-750",  
#                 "750-1000","1000-1500","1500-2000", ">2000km",  "50-200", "200-500", "500-1000")
# 
# coef_9region$dist <- factor(coef_9region$dist, levels = bin_dist, labels = label_dist)
# coef_9region <- coef_9region %>% mutate(se=se*1e2, beta=beta*1e2) 
# 
# ggplot(filter(coef_9region, model%in%c("No wind"), region!="US"), aes(x=NA, y=beta, colour=region)) +
#   geom_point(position=position_dodge(width=0.5), size=3) +
#   geom_errorbar(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), position=position_dodge(width=0.5), width=0, size=1) +
#   #geom_hline(data=filter(coef_9region, model%in%c("No wind"), region=="US"), aes(yintercept=beta), colour="purple", linetype="dashed", size=1.3) +
#   #geom_errorbar(data=filter(coef, var!="year", pval > 0.05), aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), position=position_dodge(width=0.9), colour="grey40") +
#   theme_bw() + theme(text=element_text(size=20), axis.text.x = element_text(angle=90)) +
#   labs(x="",y="ug smoke per\n 100 km2 burned area") +  #labs(x="",y="ug smoke per\n100k ton fire DM") + 
#   geom_hline(aes(yintercept=0), linetype="dashed") +
#   facet_wrap(~dist, scale="free_y", ncol=3)
# ggsave("gfedPM25_region9_pooledUS_coef_nowind.png",width=12, height = 7)
# 
# dist_list <-  c(">2000km")
# ggplot(filter(coef_9region, model%in%c("By wind (all)"), region!="US", dist %in% dist_list), aes(x=NA, y=beta, colour=region)) +
#   geom_point(position=position_dodge(width=0.5), size=3) +
#   geom_errorbar(aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), position=position_dodge(width=0.5), width=0, size=1) +
#   geom_hline(data=filter(coef_9region, model%in%c("By wind (all)"), region=="US", dist %in% dist_list), aes(yintercept=beta), colour="purple", linetype="dashed", size=1.3) +
#   #geom_errorbar(data=filter(coef, var!="year", pval > 0.05), aes(y=beta, ymin=beta-1.96*se, ymax=beta+1.96*se), position=position_dodge(width=0.9), colour="grey40") +
#   theme_bw() + theme(text=element_text(size=20), axis.text.x = element_text(angle=90)) +
#   labs(x="",y="ug smoke per\n100k ton fire DM") + geom_hline(aes(yintercept=0), linetype="dashed") +
#   facet_wrap(~dist+wind, ncol=3)
# ggsave("region9_pooledUS_coef_wind_2000km.png",width=11, height = 4)

################################    check predictions between different models (esp. pooled vs wind)    ###############################
# pred_9region <-  readRDS(paste0(result_path,"/gfedPM25_pred_grid_smoke_wind_9region_90cone.rds")) %>% mutate(type="By regions")
# pred_9region$model <- factor(pred_9region$model,  levels = c("obs","nomet","nowind","wind_main","wind_coarse","wind_all","wind_all01"),
#                              labels = c("Obs","No met, no wind","No wind","By wind","By wind (coarse)","By wind (all)","By wind (only up+down)")) 
# grid_region <- distinct(pred_9region, grid_id_10km, region)
# 
# pred_pooled <- readRDS(paste0(result_path,"/gfedPM25_pred_grid_smoke_wind_US_90cone.rds"))  %>% mutate(type="Pooled US") 
# pred_pooled$model <- factor(pred_pooled$model,  levels = c("obs","nomet","nowind","wind_main","wind_coarse","wind_all","wind_all01"),
#                             labels = c("Obs","No met, no wind","No wind","By wind","By wind (coarse)","By wind (all)","By wind (only up+down)"))
# 
# pred_pooled <- pred_pooled %>% select(-region) %>% full_join(grid_region, by="grid_id_10km")
# 
# pred_combined <- bind_rows(pred_pooled, pred_9region) %>% group_by(model, type, year, month, region) %>%
#   summarise(pred=mean(pred,na.rm=T)) %>%
  
  
