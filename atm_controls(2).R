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


setwd("C:/Users/LAPTOP/Documents")


isotopes <- read_excel("~/MUSES/1st publication_old/data/Precip&isotopes_Belvaux_2016-2023_update_02022023.xlsx", 
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



df_traj <- read.csv("~/MUSES/1st publication/new_figures/hysplits(new_Aug25).csv", header = TRUE, sep = ",", dec = ".")



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
  group_by(month, season) %>%  
  reframe(
    amount_western = round(sum(ifelse(lon < 5.94, 1, 0))/n()*100, digits=1),
    N = n()/121,
    label = paste0("N = ", N)
  ) 


ddf_Ntraj %>%
  group_by(season) %>%  
  reframe(
    N = sum(N),
    n = mean(amount_western)
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
  geom_text(data = ddf_Ntraj, aes(x = -40, y = 30, label = label), size = 2) +
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






# Define origin (same as your trajectory starting point)
origin <- c(lon = 5.94, lat = 49.51)   # <-- replace with your actual origin coordinates

# Generate circle boundaries
circle_500  <- destPoint(origin, seq(0, 360, by = 1), 500 * 1000)
circle_1500 <- destPoint(origin, seq(0, 360, by = 1), 1500 * 1000)

# Turn into data frames for ggplot
circle_500  <- as.data.frame(circle_500)
circle_1500 <- as.data.frame(circle_1500)

# Generate bearing lines (N, E, S, W, SW, etc.)
bearings <- c(45, 120, 240, 315, 200)  # match your region boundaries

lines_df <- lapply(bearings, function(b) {
  destPoint(origin, b, seq(0, 7000, by = 100) * 1000) %>%
    as.data.frame() %>%
    mutate(bearing = b)
}) %>%
  bind_rows()


region_labels <- data.frame(
  lon = c(-45, 0, 45, 20, -20),
  lat = c(35, 70, 50, 27, 27),
  label = c("W", "N", "E", "SE", "SW")
)


ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_moisture, prop_uptake*100 > 0.1), aes(color=prop_uptake*100), size=1, alpha=0.8) + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-100, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","darkblue", "purple")) +
  labs(x="", y="", color="ΔSH [%]") + theme_bw() + 
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_moisture <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_moisture, prop_uptake*100 > 0.1), aes(color=prop_uptake*100), size=1, alpha=0.9) + 
  geom_path(data = circle_500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_path(data = circle_1500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_line(data = lines_df, aes(x = lon, y = lat, group = bearing),
            color = "black", linewidth = 0.3) +
  geom_text(data = region_labels, aes(x = lon, y = lat, label = label), 
            size = 2.5, fontface= "bold") +
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","darkblue","purple", "darkred","darkred",
                                    "darkred", "darkred")) +
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
  geom_path(data = circle_500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_path(data = circle_1500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_line(data = lines_df, aes(x = lon, y = lat, group = bearing),
            color = "black", linewidth = 0.3) +
  geom_text(data = region_labels, aes(x = lon, y = lat, label = label), 
            size = 2.5, fontface= "bold") +
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  coord_sf(xlim = c(-50, 50), ylim = c(25, 75)) + 
  scale_color_gradientn(colours = c("skyblue","darkblue", "purple")) +
  labs(x="", y="", color="ΔSH [g/kg]") + theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8)) 


plot_moisture.bis








df_meteo <- read.csv("~/MUSES/1st publication_old/data/meteo_data.csv")

df_meteo$date <- as.Date(df_meteo$date)

dff_meteo <- df_meteo %>%
  mutate(temp = rollapply(temp, width = 5, FUN = mean, align = "center", fill = NA, na.rm = TRUE),
         rh = rollapply(rh, width = 5, FUN = mean, align = "center", fill = NA, na.rm = TRUE),
         gph = rollapply(gph, width = 5, FUN = mean, align = "center", fill = NA, na.rm = TRUE))

ddf <- merge(dff, dff_meteo)
ddf <- drop_na(ddf)



model <- ddf %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dO ~ temp + rh + gph + P, data = .x)),
    summary = map(model, broom::tidy),
    glance = map(model, broom::glance),
    residuals = map(model, ~ resid(.x))
  ) 


model %>%
  unnest(glance) %>%
  reframe(r_squared = round(r.squared, 2),
          sigma = round(sigma, 1),
          AIC = round(AIC, 0))


model <- ddf %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(dO ~ temp + gph + P, data = .x)),
    summary = map(model, broom::tidy),
    glance = map(model, broom::glance),
    residuals = map(model, ~ resid(.x)),
    fits = map(model, ~ fitted(.x))
  ) 


stats <- model %>%
  unnest(glance) %>%
  reframe(r_squared = round(r.squared, 2),
          sigma = round(sigma, 1),
          AIC = round(AIC, 0))

stats


summary <- model %>%
  select(summary) %>%
  unnest(summary)


ddf$fits <- unlist(model$fits)
ddf$residuals <- unlist(model$residuals)



plot_dO <- ggplot(ddf) + 
  geom_point(aes(x=date, y=dO), color="darkgrey", shape=4, size=0.5) + 
  geom_point(aes(x=date, y=fits), color="darkred", size=0.5, alpha=0.9) +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") +
  theme_bw() + labs(x="", y=expression(paste(delta^{18}, "O [\u2030]"))) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=7))




plot_dO.bis <- ggplot(ddf) + 
  geom_abline(aes(intercept=0, slope=1), linetype = "dashed") +
  geom_point(aes(x=dO, y=fits), color="darkgrey", shape=4, size=0.5) + 
  annotate("text", x = Inf, y = -Inf, hjust = 1.4, vjust = -2,
           label = paste0("R² = ", stats$r_squared, "\n", "RMSE = ", stats$sigma, " \u2030"), size = 2) +
  theme_bw() + labs(x=expression(paste("Observed",delta^{18}, "O [\u2030]")), 
                    y=expression(paste("Fitted",delta^{18}, "O [\u2030]"))) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=7)) +
  xlim(-25,5) + ylim(-25,5)



plot_dO <- ggarrange(plot_dO, plot_dO.bis,
                         labels = c("(a)", "(b)"), nrow = 1, widths = c(1.5,1),
                         font.label = list(size = 9, face = "plain"))

plot_dO








ddf2 <- dff %>%
  group_by(date, month, year) %>%
  reframe(dO = weighted.mean(dO, P, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, P, na.rm = TRUE),
          P = sum(P, na.rm=T))

ddf2 <- ddf2 %>% complete(date = seq(min(date), max(date), by = "day"))

ddf2 <- ddf2 %>%
  mutate(dO.rm = rollapply(dO*P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         d_excess.rm = rollapply(d_excess*P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         P.rm = rollapply(P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         dO.rm = dO.rm/P.rm,
         d_excess.rm = d_excess.rm/P.rm,
         dO.nm = dO - dO.rm,
         d_excess.nm = d_excess - d_excess.rm)

ddf2 <- drop_na(ddf2)

ddf2 <- merge(dff, ddf2[,c(1,7,8)])


ddf2 <- ddf2 %>%
  mutate(dO.nm = dO - dO.rm,
         d_excess.nm = d_excess - d_excess.rm)



ggplot(ddf2) + geom_point(aes(x=date, y=d_excess.nm)) + theme_bw()

ggplot(ddf2) + geom_point(aes(x=date, y=dO.nm)) + theme_bw()
ggplot(ddf2) + geom_point(aes(x=date, y=dO)) + theme_bw()









df_origin <- merge(ddf[,-c(4:6)], ddf_traj[,c(5,6,22,26:29)], by="index")



dff_origin <- df_origin %>%
  mutate(lat = floor(lat),
         lon = floor(lon)) %>%
  group_by(season) %>%
  mutate(seas_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T))


ddf_origin <- subset(dff_origin, delta_q > 0) %>%
  group_by(lat, lon, season) %>%
  reframe(uptake = sum(delta_q),
          prop_uptake = uptake/mean(seas_uptake),
          dO = weighted.mean(dO, delta_q, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, delta_q, na.rm = TRUE),
          residuals = weighted.mean(residuals, delta_q, na.rm = TRUE),
          P = sum(P, na.rm=T),
          count = n()) %>%
  ungroup()





df_origin2 <- merge(ddf2[,-c(1,3:5)], ddf_traj[,c(5,6,22,26:29)], by="index")



dff_origin2 <- df_origin2 %>%
  mutate(lat = floor(lat),
         lon = floor(lon)) %>%
  group_by(season) %>%
  mutate(seas_uptake = sum(ifelse(delta_q > 0, delta_q, 0), na.rm=T))


ddf_origin2 <- subset(dff_origin2, delta_q > 0) %>%
  group_by(lat, lon, season) %>%
  reframe(uptake = sum(delta_q),
          prop_uptake = uptake/mean(seas_uptake),
          dO = weighted.mean(dO, delta_q, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, delta_q, na.rm = TRUE),
          dO.nm = weighted.mean(dO.nm, delta_q, na.rm = TRUE),
          d_excess.nm = weighted.mean(d_excess.nm, delta_q, na.rm = TRUE),
          P = sum(P, na.rm=T),
          count = n()) %>%
  ungroup()


sf_origin2 <- st_as_sf(ddf_origin2, coords = c("lon", "lat"), crs = 4326)


region_labels <- data.frame(
  lon = c(-45, 0, 29.5, 20, -20),
  lat = c(35, 64.5, 50, 27, 27),
  label = c("W", "N", "E", "SE", "SW")
)


plot_isomap1 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_path(data = circle_500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_path(data = circle_1500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_line(data = lines_df, aes(x = lon, y = lat, group = bearing),
            color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & abs(dO.nm) >= 0.5 & abs(dO.nm) <= 3.0), 
          aes(color=dO.nm), size=1) + 
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & dO.nm > 3.0), 
          aes(), size=1.5, shape=23, fill="#C32F27", color="#C32F27") + 
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & dO.nm < -3.0), 
          aes(), size=1.5, shape=23, fill="#1F2421", color="#1F2421") + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  geom_text(data = region_labels, aes(x = lon, y = lat, label = label), 
            size = 2.5, fontface= "bold") +
  coord_sf(xlim = c(-50, 30), ylim = c(25, 65)) + 
  scale_color_gradientn(colours=c("#1F2421","#216869","#49A078","#9CC5A1","#DCE1DE",
                                  "#FFB627","#FF9505","#F26419","#C32F27")) +
  labs(x="", y="", color=expression(paste(delta^{18}, O[norm], " [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8),
        legend.position = "top") 



plot_isomap2 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = 0.3) +
  geom_path(data = circle_500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_path(data = circle_1500, aes(x = lon, y = lat),
            color = "black", linetype="dashed", linewidth = 0.3) +
  geom_line(data = lines_df, aes(x = lon, y = lat, group = bearing),
            color = "black", linewidth = 0.3) +
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & abs(d_excess.nm) >= 1.0 & abs(d_excess.nm) <= 5.0), 
          aes(color=d_excess.nm), size=1) + 
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & d_excess.nm > 5.0), 
          aes(), size=1.5, shape=23, fill="#C32F27", color="#C32F27") + 
  geom_sf(data = subset(sf_origin2, prop_uptake*100 > 0.1 & d_excess.nm < -5.0), 
          aes(), size=1.5, shape=23, fill="#1F2421", color="#1F2421") + 
  geom_sf(data = list_sf, shape = 22, fill = "red", color="black", size = 1.5) +
  geom_text(data = region_labels, aes(x = lon, y = lat, label = label), 
            size = 2.5, fontface= "bold") +
  coord_sf(xlim = c(-50, 30), ylim = c(25, 65)) + 
  scale_color_gradientn(colours=c("#1F2421","#216869","#49A078","#9CC5A1","#DCE1DE",
                                  "#FFB627","#FF9505","#F26419","#C32F27")) +
  labs(x="", y="", color=expression(paste(d-excess[norm], " [\u2030]"))) + 
  theme_bw() + facet_wrap(~season) +
  theme(strip.background = element_rect(color=NA, fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8),
        legend.position = "top")





plot_isomap <- ggarrange(plot_isomap1, plot_isomap2,
                           labels = c("(a)", "(b)"), nrow = 2,
                           font.label = list(size = 9, face = "plain"))

plot_isomap










## run only for the categorized moisture origin plots




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
    distance_km > 500 & distance_km <= 1500 & bearing_deg >= 200 & bearing_deg < 240 ~ "Southwest",
    distance_km > 500 & distance_km <= 1500 & bearing_deg >= 120 & bearing_deg < 200 ~ "Southeast",
    distance_km > 500 & distance_km <= 1500 & (bearing_deg >= 315 | bearing_deg < 45) ~ "North",
    distance_km > 500 & distance_km <= 1500 & bearing_deg >= 45 & bearing_deg < 120 ~ "East",
    distance_km > 1500 & bearing_deg >= 240 & bearing_deg < 315 ~ "Atlantic (> 1500 km)",
    distance_km > 1500 & (bearing_deg < 240 | bearing_deg >= 315) ~ "Other (> 1500 km)"))



dd_traj$region <- factor(dd_traj$region, 
                         levels = c("Other (> 1500 km)", "Atlantic (> 1500 km)",
                                    "West", "Southwest","North", "East", "Southeast",
                                    "Local (<= 500 km)"))


df_origin <- merge(df_origin, dd_traj[,c(1,13)], by="index")
ddf_gwt <- merge(ddf[,c(3,7:12,17,19,20:21)], ddf2[,c(2,18,19)], by="index")
ddf_gwt <- merge(ddf_gwt, dd_traj[,c(1:6,13)], by="index")




ddf.mm <- ddf_gwt %>%
  group_by(month, region) %>%
  reframe(P = sum(P)) %>%
  complete(month, region, fill = list(P = 0)) %>%
  group_by(month) %>%
  mutate(P_tot = sum(P),
         contr = round(P/P_tot,3),
         contr_100 = contr*100) 



ddf.mm_height <- ddf_gwt %>%
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
  geom_text(aes(x="", y=contr, label=paste0(height_i, "\nm agl")), 
            position=position_stack(vjust=0.5), size=1.5, color="white") +
  theme_minimal() + coord_polar(theta = "y") + facet_wrap(~ season) +
  labs(x="", y="", fill="Height [m agl]") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=7),
        legend.position = "none")





plot_categories <- ggarrange(plot_stacked, plot_stacked2,
                             nrow=1, widths = c(1.3,1),
                             labels = c("(a)","(b)"),
                             font.label = list(size = 8, face = "plain"))


plot_categories









ddf_gwt.mean <- ddf_gwt %>% 
  group_by(region, season) %>%
  reframe(dO_sd = round(sd(dO.nm, na.rm = TRUE),1),
          de_sd = round(sd(d_excess.nm, na.rm = TRUE),1),
          dO.nm = round(weighted.mean(dO.nm, P, na.rm = TRUE),1),
          d_excess.nm = round(weighted.mean(d_excess.nm, P, na.rm = TRUE),1),
          P = sum(P, na.rm=T))


ddf_gwt.sd <- ddf_gwt %>% 
  group_by(season) %>%
  reframe(dO_sd = round(sd(dO.nm, na.rm = TRUE),1),
          de_sd = round(sd(d_excess.nm, na.rm = TRUE),1))


ddf_gwt.mean <- arrange(ddf_gwt.mean, desc(P))
ddf_gwt.mean <- arrange(ddf_gwt.mean, season)



ddf_gwt.range <- ddf_gwt.mean %>% 
  group_by(region) %>%
  reframe(dO.nm_min = min(dO.nm, na.rm = TRUE),
          dO.nm_max = max(dO.nm, na.rm = TRUE),
          d_excess.nm_min = min(d_excess.nm, na.rm = TRUE),
          d_excess.nm_max = max(d_excess.nm, na.rm = TRUE),
          dO.nm_range = round(dO.nm_max-dO.nm_min, 1),
          dexcess.nm_range = round(d_excess.nm_max-d_excess.nm_min, 1))





plot_gwt.a <- ggplot(ddf_gwt.mean) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  geom_segment(data=ddf_gwt.range, aes(x=region, y=dO.nm_min, yend=dO.nm_max),
               linetype = "dashed") +
  geom_point(aes(x=region, y=dO.nm, fill=season, size=P), shape=21, alpha=0.9) +
  scale_fill_manual(name = "Season",
                    values = c("#03045E","#2a9d8f","#F48C06","darkgreen")) +
  scale_size_continuous(range = c(1, 12),
                        breaks = c(100, 250, 500, 1000, 1500),   
                        labels = c("0-100", "100-250", "250-500", "500-1000", "1000-1500"),
                        name = "P [mm]") +
  labs(x="Moisture origin", y=expression(paste(delta^{18}, O[norm], " [\u2030]"))) +
  theme_bw() + guides(fill = guide_legend(override.aes = list(size = 4))) + 
  theme(text = element_text(size=8)) + 
  guides(size = guide_legend(override.aes = list(shape=21))) + coord_flip()



plot_gwt.c <- ggplot(ddf_gwt.mean) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  geom_segment(data=ddf_gwt.range, aes(x=region, y=d_excess.nm_min, yend=d_excess.nm_max),
               linetype = "dashed") +
  geom_point(aes(x=region, y=d_excess.nm, fill=season, size=P), shape=21, alpha=0.9) +
  scale_fill_manual(name = "Season",
                    values = c("#03045E","#2a9d8f","#F48C06","darkgreen")) +
  scale_size_continuous(range = c(1, 12),
                        breaks = c(100, 250, 500, 1000, 1500),   
                        labels = c("0-100", "100-250", "250-500", "500-1000", "1000-1500"),
                        name = "P [mm]") +
  labs(x="Moisture origin", y=expression(paste(d-excess[norm], " [\u2030]"))) +
  theme_bw() + guides(fill = guide_legend(override.aes = list(size = 4))) + 
  theme(text = element_text(size=8)) + 
  guides(size = guide_legend(override.aes = list(shape=21))) + coord_flip()




plot_gwt <- ggarrange(plot_gwt.a, plot_gwt.c,
                             nrow=2, common.legend = TRUE, legend = "right",
                             labels = c("(a)","(b)"),
                             font.label = list(size = 8, face = "plain"))

plot_gwt





stats <- ddf_gwt %>%
  group_by(region) %>%
  reframe(x0 = median(dO, na.rm=T),
          y0 = median(fits, na.rm=T),
          median_res = round((x0 - y0)/sqrt(2),1),
          sigma = round(sqrt(sum(residuals^2)/n()),1))

labels <- stats %>%
  mutate(label = paste0("RMSE = ", round(sigma, 2), " \u2030\n",
                        "Offset = ", round(median_res, 2), " \u2030"))




plot_gwt.res <-  ggplot(ddf_gwt) +
  geom_abline(data=stats, aes(intercept=sigma/2, slope=1), linetype = "solid", color="red") +
  geom_abline(data=stats, aes(intercept=-sigma/2, slope=1), linetype = "solid", color="red") +
  geom_abline(aes(intercept=0, slope=1), linetype = "dashed") +
  geom_point(aes(x=dO, y=fits), color="#264653", shape=4, size=1) + 
  geom_point(data=stats, aes(x=x0, y=y0), color="red", shape=18, size=2) + 
  geom_text(data = labels, aes(x = 4, y = -24, label = label), size = 2, hjust = 1, vjust = 0) +
  theme_bw() + labs(x=expression(paste("Observed",delta^{18}, "O [\u2030]")), 
                    y=expression(paste("Fitted",delta^{18}, "O [\u2030]"))) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=7),
        strip.background = element_blank()) + 
  facet_wrap(~region, scales="free") + xlim(-25,5) + ylim(-25,5)


plot_gwt.res





ddf_gwt %>% 
  group_by(region) %>%
  reframe(residual.median = median(residuals, na.rm=T),
          residual.sd = sd(residuals, na.rm=T))





df_pearson <- ddf_gwt %>%
  group_by(region) %>%
  summarise(
    cor_test = list(cor.test(temp, dO, method = "pearson", use = "complete.obs")),
    pearson_rho = round(cor_test[[1]]$estimate, digits=2),
    p_value = round(cor_test[[1]]$p.value, digits=3),
    N = n()
  ) %>%
  select(region, pearson_rho, p_value, N)


df_pearson







setwd("C:/Users/LAPTOP/Documents/MUSES/1st publication")




ggsave(plot_traj, path="new_figures", filename = "hysplit_new2.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_dO, path="new_figures", filename = "mlrm_fit.png",
       device = ragg::agg_png, dpi=350,
       width = 15, height = 6, units = "cm",
       bg="white")


ggsave(plot_moisture, path="new_figures", filename = "hysplit_moisture.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_moisture.bis, path="new_figures", filename = "hysplit_moisture(abs).png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_isomap1, path="new_figures", filename = "isomap_dO(new).png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_isomap2, path="new_figures", filename = "isomap_dexcess(new).png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 16, units = "cm",
       bg="white")


ggsave(plot_gwt, path="new_figures", filename = "gwt_dO means(new).png",
       device = ragg::agg_png, dpi=350,
       width = 12, height = 12, units = "cm",
       bg="white")


ggsave(plot_categories, path="new_figures", filename = "hysplit_stacked.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 8, units = "cm",
       bg="white")


ggsave(plot_gwt.res, path="new_figures", filename = "dO_residuals.png",
       device = ragg::agg_png, dpi=350,
       width = 15, height = 15, units = "cm",
       bg="white")





write.table(summary, file='C:\\Users\\LAPTOP\\Documents\\MUSES\\1st publication\\new_figures\\summary_mlrm.csv', sep = ",", quote = FALSE)
write.table(df_pearson, file='C:\\Users\\LAPTOP\\Documents\\MUSES\\1st publication\\new_figures\\pearson_rho.csv', sep = ",", quote = FALSE)
write.table(ddf_gwt.mean, file='C:\\Users\\LAPTOP\\Documents\\MUSES\\1st publication\\new_figures\\origin_d18O.csv', sep = ",", quote = FALSE)






ggsave(cover_page, path="new_figures", filename = "cover_page.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 10, units = "cm",
       bg="white")






