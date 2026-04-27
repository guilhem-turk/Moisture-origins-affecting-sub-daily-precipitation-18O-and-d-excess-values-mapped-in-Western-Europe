library(splitr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)     # For spatial data handling
library(readxl)
library(purrr)
library(viridis)
library(patchwork)
library(zoo)
library(ggpubr)
library(lubridate)
library(ggpattern)
library(broom)

library(geosphere)  # For distance and bearing functions


setwd("C:/hysplit/metdata")


Sys.setenv(MET_DIR = "C:/hysplit/metdata") 


isotopes <- read_excel("~/MUSES/1st publication/data/Precip&isotopes_Belvaux_2016-2023_update_02022023.xlsx", 
                       col_types = c("numeric", "date", "date", 
                                     "text", "numeric", "numeric", "numeric", 
                                     "numeric", "numeric"))


isotopes <- isotopes %>% drop_na(`Sampling start`)
colnames(isotopes) <- c("index","start","end","type","P","dD","dD_dev","dO","dO_dev")


df <- isotopes %>% 
  mutate(d_excess = dD - 7.46 *dO - 5.99,
         #date = date(start), 
         date = date(start + ceiling(as.numeric(end - start)/2)), 
         time = hour(start + ceiling(as.numeric(end - start)/2)),
         year = year(date),
         month = month(date),
         duration = ceiling(as.numeric(end - start)/3600))



subset(df, P > 0 & duration <= 24) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())


subset(df, P > 1.0 & duration <= 24) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())



subset(df, P > 2.5 & duration <= 24) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())



subset(df, P > 10.0 & duration <= 24) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())




## here we select all the events of more than 1 mm in less than 24 hours

dff <- subset(df, P > 1.0 & duration <= 24)

## dff <- dff[1:5,]




## added this for date correction, don't run

df.bis <- isotopes %>% 
  mutate(d_excess = dD - 8*dO,
         date = date(start + ceiling(as.numeric(end - start)/2)), 
         time = hour(start + ceiling(as.numeric(end - start)/2)),
         year = year(date),
         month = month(date),
         duration = ceiling(as.numeric(end - start)/3600))

dff.bis <- subset(df.bis, P > 1.0 & duration <= 24)

dff.bis_helper <- merge(dff, dff.bis)
dff.bis_helper$flag <- 0

dff.bis <- merge(dff.bis, dff.bis_helper, all=T)
dff.bis$flag[is.na(dff.bis$flag)==TRUE] <- 1

ddf.bis <- subset(dff.bis, flag == 1)

#dff <- ddf.bis[,-16]


ids_to_keep <- dff.bis_helper[[1]]
df_traj.bis <- df_traj[df_traj$index %in% ids_to_keep, ]

dff_traj.bis <- rbind(df_traj.bis, df_traj.cor)

#df_traj <- dff_traj.bis


##






# Run HYSPLIT safely and bind results


results_500 <- pmap_dfr(
  list(seq_along(dff$index), dff$index, dff$date, dff$time),
  function(i, index, date, time) {
    total <- nrow(dff)
    percent_done <- round((i / total) * 100)
    
    message(sprintf(
      "Running trajectory %d of %d (%d%%) - Index: %s | Date: %s at %sh",
      i, total, percent_done, index, date, time
    ))
    
    tryCatch({
      traj <- hysplit_trajectory(
        lat = 49.62,
        lon = 6.22,
        height = 500,
        duration = 120,
        days = date,
        daily_hours = time,
        met_type = "reanalysis",
        extended_met = TRUE,
        direction = "backward"
      )
      traj$index <- index
      traj
    }, error = function(e) {
      warning("Failed on ", index, " - ", date, " at ", time, ": ", e$message)
      NULL
    })
  }
)


results_1000 <- pmap_dfr(
  list(seq_along(dff$index), dff$index, dff$date, dff$time),
  function(i, index, date, time) {
    total <- nrow(dff)
    percent_done <- round((i / total) * 100)
    
    message(sprintf(
      "Running trajectory %d of %d (%d%%) - Index: %s | Date: %s at %sh",
      i, total, percent_done, index, date, time
    ))
    
    tryCatch({
      traj <- hysplit_trajectory(
        lat = 49.62,
        lon = 6.22,
        height = 1000,
        duration = 120,
        days = date,
        daily_hours = time,
        met_type = "reanalysis",
        extended_met = TRUE,
        direction = "backward"
      )
      traj$index <- index
      traj
    }, error = function(e) {
      warning("Failed on ", index, " - ", date, " at ", time, ": ", e$message)
      NULL
    })
  }
)


results_1500 <- pmap_dfr(
  list(seq_along(dff$index), dff$index, dff$date, dff$time),
  function(i, index, date, time) {
    total <- nrow(dff)
    percent_done <- round((i / total) * 100)
    
    message(sprintf(
      "Running trajectory %d of %d (%d%%) - Index: %s | Date: %s at %sh",
      i, total, percent_done, index, date, time
    ))
    
    tryCatch({
      traj <- hysplit_trajectory(
        lat = 49.62,
        lon = 6.22,
        height = 1500,
        duration = 120,
        days = date,
        daily_hours = time,
        met_type = "reanalysis",
        extended_met = TRUE,
        direction = "backward"
      )
      traj$index <- index
      traj
    }, error = function(e) {
      warning("Failed on ", index, " - ", date, " at ", time, ": ", e$message)
      NULL
    })
  }
)



results_2000 <- pmap_dfr(
  list(seq_along(dff$index), dff$index, dff$date, dff$time),
  function(i, index, date, time) {
    total <- nrow(dff)
    percent_done <- round((i / total) * 100)
    
    message(sprintf(
      "Running trajectory %d of %d (%d%%) - Index: %s | Date: %s at %sh",
      i, total, percent_done, index, date, time
    ))
    
    tryCatch({
      traj <- hysplit_trajectory(
        lat = 49.62,
        lon = 6.22,
        height = 2000,
        duration = 120,
        days = date,
        daily_hours = time,
        met_type = "reanalysis",
        extended_met = TRUE,
        direction = "backward"
      )
      traj$index <- index
      traj
    }, error = function(e) {
      warning("Failed on ", index, " - ", date, " at ", time, ": ", e$message)
      NULL
    })
  }
)




write.table(results_500, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits_500.csv', sep = ",", quote = FALSE)
write.table(results_1000, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits_1000.csv', sep = ",", quote = FALSE)
write.table(results_1500, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits_1500.csv', sep = ",", quote = FALSE)
write.table(results_2000, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits_2000.csv', sep = ",", quote = FALSE)


df_traj <- rbind(results_500, results_1000) 
df_traj <- rbind(df_traj, results_1500) 
df_traj <- rbind(df_traj, results_2000)  




write.table(df_traj, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits(all_layers).csv', sep = ",", quote = FALSE)


##


## df_traj <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/new_figures/hysplits(all_layers).csv", header = TRUE, sep = ",", dec = ".")

df_traj <- read.csv("C:/Users/turk/Documents/MUSES/Isotope_model/new_figures/hysplits(new_Aug25).csv", header = TRUE, sep = ",", dec = ".")





##


## once the trajectories have been calculated and stored, I recommend using atm_controls(2).R
## it contains the most recent versions of the figures used in the manuscript



##























df_traj <- df_traj  %>% 
  mutate(date = date(traj_dt_i), 
         year = year(date),
         month = month(date))


get_season <- function(date) {
  month = month(date)
  if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else if (month %in% c(9, 10, 11)) {
    return("Autumn")
  } else {
    return("Winter")
  }
}

df_traj$season <- sapply(df_traj$date, get_season)

df_traj$season <- factor(df_traj$season, levels = c("Winter", "Spring", "Summer", "Autumn"))


df_Ntraj <- df_traj %>%
  group_by(index) %>%  
  reframe(
    traj_end = min(hour_along)
  ) 


df_traj <- df_traj %>%
  group_by(index, height_i) %>%  
  mutate(
    traj_end = min(hour_along)
  ) 


dff_traj <- subset(df_traj, traj_end <= -24) %>%
  group_by(index, height_i) %>%  
  arrange(hour_along, .by_group = TRUE) %>%
  mutate(
    delta_q = sp_humidity - lag(sp_humidity),  
    delta_q6d = rollapply(delta_q, width = 6, FUN = sum, align = "right", fill = NA, na.rm = TRUE),
    uptake = ifelse(delta_q6d >= 0.5, 1, 0),
    tot_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T),
    prop_uptake = ifelse(delta_q > 0, delta_q/tot_uptake, 0)
  ) %>%
  ungroup()



ddf_traj <- dff_traj %>%
  group_by(index) %>%
  filter(tot_uptake == max(tot_uptake))


ddf_Ntraj <- ddf_traj %>%
  group_by(index) %>%  
  reframe(
    traj_end = min(hour_along)
  ) 


ddf_traj$month_name <- months(ddf_traj$date, abbreviate=T)

ddf_traj$month_name <- factor(ddf_traj$month_name, 
                            levels=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug", "Sep","Oct","Nov"))



sf_traj <- st_as_sf(drop_na(ddf_traj), coords = c("lon", "lat"), crs = 4326)

world_map <- map_data("world")

list_sf <- st_sf(
  geometry = st_sfc(st_point(c(5.94, 49.51)), crs = 4326)
)








plot_traj <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = sf_traj, aes(), color="grey", size=0.2, shape=16, alpha=0.7) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  labs(x="", y="") + theme_bw() + facet_wrap(~month_name) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_traj





cover_page <- ggplot(sf_traj) +
  geom_sf(aes(color=hour_along), size=0.4, shape=16, alpha=0.7) + 
  scale_color_viridis(option="magma") +
  coord_sf(xlim = c(-90, 20), ylim = c(35, 65)) + 
  labs(x="", y="") + theme_minimal() + 
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=5),
        legend.position = "none") 


cover_page








ddf_traj %>%
  group_by(season) %>%
  reframe(seas_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T))


dd_traj <- ddf_traj %>%
  mutate(lat = floor(lat),
         lon = floor(lon)) %>%
  group_by(season) %>%
  mutate(seas_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T))


dd_traj <- dd_traj %>%
  group_by(lat, lon, season) %>%
  reframe(uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T),
          prop_uptake = uptake/mean(seas_uptake),
          count = n()) %>%
  ungroup()





sf_moisture <- st_as_sf(dd_traj, coords = c("lon", "lat"), crs = 4326)





plot_moisture <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_moisture, prop_uptake*100 > 0.1), aes(color=prop_uptake*100), size=1, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","darkblue", "purple")) +
  labs(x="", y="", color="ΔSH [%]") + theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_moisture




plot_moisture.bis <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_moisture, prop_uptake*100 > 0.1), aes(color=uptake), size=1, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","darkblue", "purple")) +
  labs(x="", y="", color="ΔSH [g/kg]") + theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_moisture.bis












df_origin <- merge(dff[,-c(2:4)], ddf_traj[,c(5,6,22,26:29)], by="index")


dff_origin <- df_origin %>%
  mutate(lat = floor(lat/2)*2,
         lon = floor(lon/2)*2) %>%
  group_by(season) %>%
  mutate(seas_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T))


ddf_origin <- subset(dff_origin, delta_q > 0) %>%
  group_by(lat, lon, season) %>%
  reframe(uptake = sum(delta_q),
          prop_uptake = uptake/mean(seas_uptake),
          dO = weighted.mean(dO, delta_q, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, delta_q, na.rm = TRUE),
          P = sum(P, na.rm=T),
          count = n()) %>%
  ungroup()





sf_origin <- st_as_sf(ddf_origin, coords = c("lon", "lat"), crs = 4326)





plot_origin <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2), aes(color=dO), size=1, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color=expression(paste(delta^{18}, "O [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_origin





plot_origin2 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2), aes(color=d_excess), size=1, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color="lc-excess [\u2030]") + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_origin2










dd_traj <- subset(ddf_traj, delta_q > 0) %>%
  group_by(index, date, year, month, season, height_i) %>%
  reframe(lat = weighted.mean(lat, delta_q),
          lon = weighted.mean(lon, delta_q),
          uptake = sum(delta_q),
          count = n()) %>%
  ungroup()




receptor <- c(lon = 6.22, lat = 49.62)  # Your sampling site


dd_traj <- dd_traj %>%
  mutate(
    distance_km = distHaversine(cbind(lon, lat), receptor) / 1000,
    bearing_deg = (bearing(receptor, cbind(lon, lat)) + 360) %% 360
  )



dd_traj <- dd_traj %>%
  mutate(region = case_when(
      distance_km <= 500 ~ "Local (<= 500 km)",
      distance_km > 500 & distance_km <= 1500 & bearing_deg >= 225 & bearing_deg < 315 ~ "West",
      distance_km > 500 & distance_km <= 1500 & bearing_deg >= 200 & bearing_deg < 225 ~ "Southwest",
      distance_km > 500 & distance_km <= 1500 & bearing_deg >= 135 & bearing_deg < 200 ~ "South",
      distance_km > 500 & distance_km <= 1500 & (bearing_deg >= 315 | bearing_deg < 45) ~ "North",
      distance_km > 500 & distance_km <= 1500 & bearing_deg >= 45 & bearing_deg < 135 ~ "East",
      distance_km > 1500 & bearing_deg >= 225 & bearing_deg < 315 ~ "West (> 1500 km)",
      distance_km > 1500 & (bearing_deg < 225 | bearing_deg >= 315) ~ "Other (> 1500 km)"))



dd_traj$region <- factor(dd_traj$region, 
                         levels = c("Other (> 1500 km)", "West (> 1500 km)",
                                    "West", "Southwest","North", "East", "South",
                                    "Local (<= 500 km)"))


df_origin <- merge(df_origin, dd_traj[,c(1,13)], by="index")
ddf <- merge(dff[,c(1,5:10)], dd_traj[,c(1:6,13)], by="index")




ddf.mm <- ddf %>%
  group_by(month, region) %>%
  reframe(P = sum(P)) %>%
  complete(month, region, fill = list(P = 0)) %>%
  group_by(month) %>%
  mutate(P_tot = sum(P),
         contr = P/P_tot) 



ddf.mm_height <- ddf %>%
  group_by(season, height_i) %>%
  reframe(P = sum(P)) %>%
  complete(season, height_i, fill = list(P = 0)) %>%
  group_by(season) %>%
  mutate(P_tot = sum(P),
         contr = P/P_tot) 




plot_stacked <- ggplot(ddf.mm) + 
  geom_area(aes(x=month, y=contr, fill=region),
                    alpha=0.8, position="stack", color="white", lwd=0.2) + 
  theme_minimal() + scale_fill_manual(values = c("#264653","darkgrey","#03045E","#2a9d8f",
                                                 "skyblue","#e9c46a","#f4a261","lightgrey")) +
  scale_x_continuous(breaks = 1:12, labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  labs(x="Month", y="Fraction of total rain amount [-]", fill="Moisture origin") +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=7))




plot_stacked2 <- ggplot(ddf.mm_height) + 
  geom_col(aes(x="", y=contr, fill=height_i), alpha=0.9, color="white", lwd=0.2) + 
  theme_minimal() + coord_polar(theta = "y") + facet_wrap(~ season) +
  labs(x="", y="", fill="Height [m asl]") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=7))





plot_origin3 <- ggarrange(plot_stacked, plot_stacked2,
                         nrow=1, widths = c(1.3,1),
                         labels = c("(a)","(b)"),
                         font.label = list(size = 8, face = "plain"))
  

plot_origin3







ggplot(ddf) +
  geom_boxplot(aes(x=region, y=dO, fill=region), alpha=0.8) +
  scale_fill_manual(values = c("#264653","darkgrey","#03045E","#2a9d8f",
                               "skyblue","#e9c46a","#f4a261","lightgrey")) +
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), fill="Moisture origin") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=8))







df_rainout <- subset(ddf_traj, hour_along < -5) %>%
  group_by(index) %>%
  reframe(rainout = sum(rainfall),
          count = n()) %>%
  ungroup()



run_model_for_window <- function(window_size) {
  
  # Filter and summarise for the given hour_along window
  df_remote <- subset(ddf_traj, hour_along >= window_size) %>%
    group_by(index, date, year, month, season) %>%
    reframe(
      temp_rm = mean(air_temp),
      rh_rm = mean(rh),
      rainout = sum(rainfall),
      count = n()
    ) %>%
    ungroup()
  
  # Merge everything
  ddf <- merge(dff[,c(1,5:10)], dd_traj[,c(1:6,13)], by="index")
  ddf <- merge(ddf, df_remote[,c(1,6,7,8)], by="index")
  
  # Fit the models
  model <- ddf %>%
    group_by(region) %>%
    nest() %>%
    mutate(
      model = map(data, ~ lm(dO ~ temp_rm + rh_rm + rainout, data = .x)),
      glance = map(model, broom::glance)
    )
  
  # Return metrics
  model %>%
    unnest(glance) %>%
    group_by(region) %>%
    reframe(
      window = window_size,
      r_squared = round(r.squared, 2),
      sigma = round(sigma, 1),
      AIC = round(AIC, 0)
    )
}


results <- map_dfr(c(-3,-6, -12, -18, -24, -36, -48, -72, -96), run_model_for_window)





plot_mlrm <- ggplot(results) +
  geom_line(aes(x=window, y=r_squared), linetype="dashed") +
  geom_point(aes(x=window, y=r_squared, color=sigma)) +
  scale_color_gradient(low = "darkblue", high = "red") +
  theme_bw() + facet_wrap(~region, scales="free") + ylim(0,1) +
  labs(x="Trajectory duration [h]", y="R² [-]", 
       color=expression(paste("RMSE ", delta^{18}, "O [\u2030]"))) +
  theme(strip.background = element_blank(),
         text = element_text(size=8))


plot_mlrm




ddf <- merge(dff[,c(1,5:10)], dd_traj[,c(1:6,13)], by="index")


df_meteo <- read.csv("~/MUSES/1st publication/data/meteo_data.csv")

df_meteo$date <- as.Date(df_meteo$date)

dff_meteo <- merge(df_meteo, dff[,c(1,11,12)])

ddf <- merge(ddf, dff_meteo[,c(4:7)], by = "index")





model <- ddf %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dO ~ temp + rh + gph + P, data = .x)),
    summary = map(model, broom::tidy),
    glance = map(model, broom::glance)
  ) 

model %>%
  unnest(glance) %>%
  reframe(r_squared = round(r.squared, 2),
          sigma = round(sigma, 1),
          AIC = round(AIC, 0))





model <- ddf %>%
  group_by(region) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dO ~ temp + rh + gph + P, data = .x)),
    summary = map(model, broom::tidy),
    glance = map(model, broom::glance)
  ) 

df_model <- model %>%
  unnest(glance) %>%
  group_by(region) %>%
  reframe(r_squared = round(r.squared, 2),
          sigma = round(sigma, 1),
          AIC = round(AIC, 0))

df_model




ggplot(results) +
  geom_line(aes(x=window, y=r_squared), linetype="dashed") +
  geom_point(aes(x=window, y=r_squared, color=sigma)) +
  geom_point(data=df_model, aes(x=0, y=r_squared, color=sigma), shape=18, size=2) +
  scale_color_gradient(low = "darkblue", high = "red") +
  theme_bw() + facet_wrap(~region, scales="free") + ylim(0,1) +
  labs(x="Trajectory duration [h]", y="R² [-]", 
       color=expression(paste("RMSE ", delta^{18}, "O [\u2030]"))) +
  theme(strip.background = element_blank(),
        text = element_text(size=8))






ggplot(ddf) +
  geom_point(aes(x=temp, y=dO)) +
  scale_color_viridis(discrete=T) +
  theme_bw()


ggplot(ddf) +
  geom_point(aes(x=rh, y=dO)) +
  scale_color_viridis(discrete=T) +
  theme_bw() 


ggplot(ddf) +
  geom_point(aes(x=gph, y=dO)) +
  scale_color_viridis(discrete=T) +
  theme_bw()










setwd("C:/Users/turk/Documents/MUSES/Isotope_model")



ggsave(plot_traj, path="new_figures", filename = "hysplit_new2.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_moisture, path="new_figures", filename = "hysplit_moisture.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_moisture.bis, path="new_figures", filename = "hysplit_moisture(abs).png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_origin, path="new_figures", filename = "hysplit_origin.png",
       device = ragg::agg_png, dpi=350,
       width = 22, height = 14, units = "cm",
       bg="white")


ggsave(plot_origin2, path="new_figures", filename = "hysplit_origin_dO.png",
       device = ragg::agg_png, dpi=350,
       width = 22, height = 14, units = "cm",
       bg="white")


ggsave(plot_origin3, path="new_figures", filename = "hysplit_stacked.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 8, units = "cm",
       bg="white")


ggsave(plot_mlrm, path="new_figures", filename = "hysplit_mlrm.png",
       device = ragg::agg_png, dpi=350,
       width = 15, height = 13, units = "cm",
       bg="white")






ggsave(cover_page, path="new_figures", filename = "cover_page.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 10, units = "cm",
       bg="white")






#write.table(df_traj, file='C:\\Users\\turk\\Documents\\MUSES\\Isotope_model\\new_figures\\hysplits(P2,5-120h).csv', sep = ",", quote = FALSE)










































## Archive










ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = sf_traj, aes(), color="grey", size=0.2, shape=16, alpha=0.4) + 
  geom_sf(data = subset(sf_traj, delta_q > 0), aes(color=delta_q), size=0.2, shape=16, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 24, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","purple","darkblue")) +
  labs(x="", y="", color="ΔSH [g/kg]") + theme_bw() + 
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 





ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = sf_traj, aes(), color="grey", size=0.2, shape=16, alpha=0.3) + 
  geom_sf(data = subset(sf_traj, delta_q > 0), aes(color=delta_q), size=0.2, shape=16, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 24, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-20, 20), ylim = c(35, 60)) + 
  scale_color_gradientn(colours = c("skyblue","purple","darkblue")) +
  labs(x="", y="", color="ΔSH [g/kg]") + theme_bw() + 
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 

















world <- ne_countries(scale = "medium", returnclass = "sf")


ocean <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")


bbox <- st_bbox(c(xmin = -50, xmax = 50, ymin = 20, ymax = 75), crs = st_crs(world))
world <- st_crop(world, bbox)
ocean <- st_crop(ocean, bbox)


world_centroids <- st_centroid(world)


coords <- st_coordinates(world_centroids)
world$longitude <- coords[,1]
world$lat <- coords[,2]



world <- world %>%
  mutate(region = case_when(
    lat < 45 ~ "SE-EA",
    lat >= 45 & longitude <= 15 ~ "WE",
    lat >= 45 & longitude > 15 ~ "EE"
  ))





bbox_poly <- st_as_sfc(st_bbox(c(xmin = -50, xmax = 50, ymin = 20, ymax = 75), crs = st_crs(ocean)))


ocean_clipped <- st_intersection(ocean, bbox_poly)



ocean_centroids <- st_centroid(ocean)

coords_ocean <- st_coordinates(ocean_centroids)
ocean$longitude <- coords_ocean[,1]
ocean$lat <- coords_ocean[,2]

ocean <- ocean %>%
  mutate(region2 = case_when(
    lat < 40 & longitude <= 5 ~ "SA",
    lat < 45 & longitude > 5 ~ "ME",
    lat >= 40 & lat < 55 & longitude < 0 ~ "AT",
    lat >= 55  ~ "NA-NS"
  ))



ggplot() +
  #geom_sf(data = world, aes(fill = region), color = "black") +
  geom_sf(data = ocean, aes(fill = region2), color = "blue", alpha = 0.4) +
  #coord_sf(xlim = c(-50, 50), ylim = c(20, 75)) +
  theme_minimal() +
  labs(title = "Regional Split: Land & Ocean")











plot_a <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "winter"), 
          aes(color=d_excess), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color="d-excess [\u2030]") + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_b <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "spring"), 
          aes(color=d_excess), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color="d-excess [\u2030]") + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_c <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "summer"), 
          aes(color=d_excess), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color="d-excess [\u2030]") + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_d <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "autumn"), 
          aes(color=d_excess), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color="d-excess [\u2030]") + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 




plot_origin <- (plot_a | plot_b) / (plot_c | plot_d) 

plot_origin










plot_a <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "winter"), 
          aes(color=dO), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color=expression(paste(delta^{18}, "O [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_b <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "spring"), 
          aes(color=dO), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color=expression(paste(delta^{18}, "O [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_c <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "summer"), 
          aes(color=dO), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color=expression(paste(delta^{18}, "O [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 



plot_d <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin, prop_uptake*100 > 0.2 & season == "autumn"), 
          aes(color=dO), size=0.7, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("yellow","red")) +
  labs(x="", y="", color=expression(paste(delta^{18}, "O [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 




plot_origin2 <- (plot_a | plot_b) / (plot_c | plot_d) 

plot_origin2






dd_traj <- dd_traj %>%
  mutate(region = case_when(
    bearing_deg >= 225 & bearing_deg < 315 ~ "West",
    bearing_deg >= 200 & bearing_deg < 225 ~ "Southwest",
    bearing_deg >= 135 & bearing_deg < 200 ~ "South",
    bearing_deg >= 315 | bearing_deg < 45 ~ "North",
    bearing_deg >= 45 & bearing_deg < 135 ~ "East"),
    distance = case_when(
      distance_km <= 500 ~ "Local",
      distance_km > 500 & distance_km <= 1500 ~ "Near",
      distance_km > 1500 ~ "Far")
  )

dd_traj$region <- factor(dd_traj$region, levels = c("West", "Southwest", "North", "East", "South"))
dd_traj$distance <- factor(dd_traj$distance, levels = c("Local", "Near", "Far"))





ddf.mm_dist <- ddf %>%
  group_by(month, distance) %>%
  reframe(P = sum(P)) %>%
  complete(month, distance, fill = list(P = 0)) %>%
  group_by(month) %>%
  mutate(P_tot = sum(P),
         contr = P/P_tot) 





plot_stacked2 <- ggplot(ddf.mm_dist) + 
  geom_area(aes(x=month, y=contr, fill=distance), alpha=0.9, position="stack", color="white", lwd=0.2) + 
  theme_minimal() +
  scale_x_continuous(breaks = 1:12, labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  labs(x="Month", y="Fraction of total rain amount [-]", fill="Trajectory") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size=7))

plot_stacked2






ddf_monthly <- ddf %>%
  group_by(year, month, season) %>%
  reframe(dO = weighted.mean(dO, P, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, P, na.rm = TRUE),
          temp_rm = weighted.mean(temp_rm, P, na.rm = TRUE),
          rh_rm = weighted.mean(rh_rm, P, na.rm = TRUE),
          rainout = weighted.mean(rainout, P, na.rm = TRUE),
          P = sum(P, na.rm=T))












df_remote <- subset(ddf_traj, hour_along >= -12) %>%
  group_by(index, date, year, month, season) %>%
  reframe(temp_rm = mean(air_temp),
          rh_rm = mean(rh),
          count = n()) %>%
  ungroup()




ddf <- merge(dff[,c(1,5:10)], dd_traj[,c(1:6,13)], by="index")

ddf <- merge(ddf, df_remote[,c(1,6,7)], by="index")
ddf <- merge(ddf, df_rainout[,c(1,2)], by="index")



model <- ddf %>%
  group_by(region) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dO ~ temp_rm + rh_rm + rainout, data = .x)),
    summary = map(model, broom::tidy),
    glance = map(model, broom::glance)
  ) 

model %>%
  unnest(glance) %>%
  group_by(region) %>%
  reframe(r_squared = round(r.squared, 2),
          sigma = round(sigma, 1),
          AIC = round(AIC, 0))



ggplot(ddf) +
  geom_point(aes(x=rainout, y=dO)) +
  theme_bw() 


ggplot(ddf) +
  geom_point(aes(x=temp_rm, y=dO)) +
  scale_color_viridis(discrete=T) +
  theme_bw()


ggplot(ddf) +
  geom_point(aes(x=rh_rm, y=dO)) +
  scale_color_viridis(discrete=T) +
  theme_bw() 


