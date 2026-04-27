library(readxl)
library(ggpubr)
library(scales)
library(dplyr)
library(tidyr)
library(stats)
library(gridExtra)
library(purrr)
library(viridis)
library(imputeTS)
library(lubridate)
library(data.table)
library(lmodel2)
library(zoo)


setwd("C:/Users/LAPTOP/Documents")


df_meteo <- read.csv("~/MUSES/1st publication_old/data/meteo_data.csv")


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


df$date <- as.Date(df$date)





dff <- subset(df, P > 1.0 & duration <= 24) %>%
  group_by(date, month, year) %>%
  reframe(dO = weighted.mean(dO, P, na.rm = TRUE),
          d_excess = weighted.mean(d_excess, P, na.rm = TRUE),
          P = sum(P, na.rm=T))



dff <- dff %>% complete(date = seq(min(date), max(date), by = "day"))


dff <- dff %>%
  mutate(dO = rollapply(dO*P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         d_excess = rollapply(d_excess*P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         P = rollapply(P, width = 30, FUN = sum, align = "center", fill = NA, na.rm = TRUE),
         dO = dO/P,
         d_excess = d_excess/P)






plot_a <- ggplot() + 
  geom_point(data=df, aes(x=date, y=dO), color="darkgrey", shape=4, size=0.6) + 
  geom_line(data=dff, aes(x=date, y=dO), color="darkred", lwd = 0.4) +  
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + theme_bw() +
  labs(x="", y=expression(paste(delta^{18}, "O [\u2030]")), linetype="") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=9))



plot_a.bis <- ggplot() + 
  geom_point(data=df, aes(x=date, y=d_excess), color="darkgrey", shape=4, size=0.6) +  
  geom_line(data=dff, aes(x=date, y=d_excess), color="darkred", lwd = 0.4) +  
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + theme_bw() +
  labs(x="", y="d-excess [\u2030]", linetype="") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=9))



dff_meteo <- subset(df_meteo, date > "2016-12-31") %>%
  mutate(year=year(date), month=month(date)) %>%
  group_by(date) %>%
  reframe(temp.sd = sd(temp, na.rm=T),
          temp = median(temp, na.rm=T),
          rh.sd = sd(rh, na.rm=T),
          rh = median(rh, na.rm=T),
          gph.sd = sd(gph, na.rm=T),
          gph = median(gph, na.rm=T)) 



dff_meteo$date <- as.Date(dff_meteo$date)



plot_b <- ggplot(dff_meteo) +
  #geom_ribbon(aes(x=date, y=temp-40, ymin=temp-40-temp.sd, ymax=temp-40+temp.sd), fill="grey", alpha=0.5) +
  geom_col(data=df, aes(x=date, y=-(P)/2), fill="skyblue", width=5) +
  geom_line(aes(x=date, y=temp-40), lwd = 0.2) +  
  scale_y_continuous(labels = function(y) y + 40, 
                     sec.axis = sec_axis(~.*(-2), name="Precipitation [mm]")) + 
  theme_bw() + labs(x="", y="Temperature [°C]") +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=9))






plot_c <- ggplot(dff_meteo) + 
  geom_line(aes(x=date, y=gph-1000), color="darkgrey", lwd = 0.2) +
  geom_line(aes(x=date, y=rh), color="darkblue", lwd = 0.2) +  
  scale_y_continuous(sec.axis = sec_axis(~.+1000, name="Geopotential height [hPa]")) + 
  theme_bw() + labs(x="", y="Relative humidity [%]") +
  scale_x_date(date_breaks="1 year", date_labels="%b %Y") + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_line(linetype = "dashed", color="black"),
        panel.grid.minor = element_blank(),
        text = element_text(size=9))





plot_comb <- ggarrange(plot_a, plot_a.bis, plot_b,  plot_c,
                       align="hv", common.legend = T,  ncol = 1, nrow = 4,
                       labels = c("(a)", "(b)",  "(c)", "(d)"), font.label = list(size = 9, face = "plain"),
                       legend = "none")


plot_comb





##



# Fit RMA model (doesn't support weights directly)
rma_model <- lmodel2(dD ~ dO, data = drop_na(df_iso), range.y = "interval", range.x = "interval", nperm = 0)

rma_results <- rma_model$regression.result

# Filter for RMA method
label <- rma_results %>%
  filter(Method == "RMA") %>%
  summarise(
    intercept = Intercept,
    slope     = Slope,
    r_squared = rma_model$rsquare,
    label = paste0(
      "d²H = ", round(Intercept, 2), " + ", round(Slope, 2), " * d¹⁸O",
      "\nR² = ", round(rma_model$rsquare, 2)
    )
  )




plot_lmwl.a <- ggplot(subset(df_iso, P >= 1.0)) +
  geom_abline(aes(intercept=10, slope=8), color="black", linetype="solid", lwd=0.3) +
  geom_abline(data=label, aes(intercept=intercept, slope=slope), color="black", linetype="dashed", lwd=0.3) +
  geom_point(data=subset(df_iso, P < 1.0), aes(x=dO, y=dD), 
             color="grey", alpha=0.7, shape=4, size=0.7) + 
  geom_point(aes(x=dO, y=dD), color="black",
             alpha=0.7, shape=4, size=0.7) + 
  geom_text(data=label, aes(x=-15, y=10, label = label), size=2) +
  labs(x=expression(paste(delta^{18}, "O [\u2030]")), y=expression(paste(delta^{2}, "H [\u2030]"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=8))






rma_models <- drop_na(df_iso) %>%
  group_by(season) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lmodel2(dD ~ dO, data = ., range.y = "interval", range.x = "interval", nperm = 0)),
    results = map(model, ~ .x$regression.result),
    rsq = map_dbl(model, ~ .x$rsquare)
  ) %>%
  select(season, results, rsq) %>%
  unnest(results) %>%
  filter(Method == "RMA") %>%
  transmute(
    season,
    intercept = Intercept,
    slope = Slope,
    r_squared = rsq,
    label = paste0(
      "d²H = ", round(Intercept, 2), " + ", round(Slope, 2), " * d¹⁸O",
      "\nR² = ", round(r_squared, 2)
    )
  )


label.pos <- drop_na(df_iso) %>%
  group_by(season) %>%
  reframe(dO = min(dO) * 0.6,
          dD = max(dD) * 0.4)


label.seas <- merge(label.pos, rma_models)



plot_lmwl.b <- ggplot(subset(df_iso, P >= 1.0)) +
  geom_abline(aes(intercept=10, slope=8), color="black", linetype="solid", lwd=0.3) +
  geom_abline(data=label.seas, aes(intercept=intercept, slope=slope), color="black", linetype="dashed", lwd=0.3) +
  geom_point(data=subset(df_iso, P < 1.0), aes(x=dO, y=dD), 
             color="grey", alpha=0.7, shape=4, size=0.4) + 
  geom_point(aes(x=dO, y=dD), color="black",
             alpha=0.7, shape=4, size=0.4) + 
  geom_text(data=label.seas, aes(x=dO, y=-10, label = label), size=2) +
  labs(x=expression(paste(delta^{18}, "O [\u2030]")), y=expression(paste(delta^{2}, "H [\u2030]"))) +
  theme_bw() + facet_wrap(~season, scales="free") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size=8)) 





plot_lmwl.comb <- ggarrange(plot_lmwl.a, plot_lmwl.b,  
                            labels = c("(a)", "(b)"), ncol=2, 
                            font.label = list(size = 8, face = "plain"))


plot_lmwl.comb





##





isotopes <- isotopes %>% drop_na(`Sampling start`)
colnames(isotopes) <- c("index","start","end","type","P","dD","dD_dev","dO","dO_dev")


isotopes <- isotopes %>% 
  mutate(d_excess = dD - 8*dO,
         date = date(start), 
         time = hour(start),
         year = year(date),
         month = month(date))

isotopes$duration <- ceiling(as.numeric(isotopes$end - isotopes$start)/3600)



subset(isotopes, P>0) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(weighted.mean(dO, P, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T),1),
          precip = sum(P),
          N = n())

df %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(weighted.mean(dO, P, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T),1),
          precip = sum(P),
          N = n())


df %>%
  group_by(season) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(weighted.mean(dO, P, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T),1),
          precip = sum(P),
          N = n())




df %>%
  group_by(region_id) %>%
  reframe(dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())

df %>%
  reframe(dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(median(dO, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(median(d_excess, na.rm=T),1),
          precip = sum(P),
          N = n())


df.mm <- df %>%
  group_by(year, month) %>%
  reframe(dO.sd = round(sd(dO, na.rm=T),1),
          dO.iqr = round(IQR(dO, na.rm=T),1),
          dO = round(weighted.mean(dO, P, na.rm=T),1),
          d_excess.sd = round(sd(d_excess, na.rm=T),1),
          d_excess.iqr = round(IQR(d_excess, na.rm=T),1),
          d_excess = round(weighted.mean(d_excess, P, na.rm=T),1),
          precip = sum(P),
          N = n())



df.mm %>%
  reframe(dO.sd = round(mean(dO.sd, na.rm=T),1),
          d_excess.sd = round(mean(d_excess.sd, na.rm=T),1),
          N = n())








setwd("C:/Users/LAPTOP/Documents/MUSES/1st publication")



ggsave(plot_comb, path="new_figures", filename = "variables.png",
       device = ragg::agg_png, dpi=350,
       width = 16, height = 22, units = "cm",
       bg="white")


ggsave(plot_lmwl.comb, path="new_figures", filename = "lmwl.png",
       device = ragg::agg_png, dpi=350,
       width = 18, height = 10, units = "cm",
       bg="white")
