#' R Code accompanying research article
#' Grundlehner et al. (2026) 
#' TITLE
#' JOURNAl
#' DOI:
#'  -- (1)  Cluster analysis using Ripleys K
#' Script by: Anne Grundlehner (anne.grundlehner@utas.edu.au), 
#' Affiliation: Institute of Marine & Antarctic studies, Hobart, Tasmania, Australia
#' 
#' **** ANY USE OF THIS CODE SHOULD REFER TO GRUNDLEHNER ET AL (2026) ****


rm(list=ls())

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(sf)
library(raster)
library(sp)
library(spatstat)
library(mgcv)

#' Set your working directory
my_wd = "~/Documents/PhD/DATA/CLEAN"

#### 0. Import data ####

setwd(my_wd) 
fg = read_xlsx("fg_grid0.1_AggrAreaAttributed.xlsx")

df0 = fg%>%
  dplyr::select(LonBin, LatBin, Year, total)%>%
  mutate(Count = total)

#' Make table to store values
my_list = df0%>%
  dplyr::select(Year)%>%
  distinct()%>%
  mutate(Ripley = NA)


#' Estimation of Ripley's K function for inhomogeneous spatial patterns
#' This is performed for each year in my_list seperately
for(i in 1:nrow(my_list)){
  
  year=my_list[i,]$Year
  df <- as.data.frame(lapply(df0, rep, df0$Count))%>% # one row for each count
    filter(Year %in% c(year))
  mygeom = st_geometry(st_as_sf(x=df, coords = c("LonBin", "LatBin")))
  mygeom$geometry = st_as_sf(x=df,  coords = c("LonBin", "LatBin"), crs=crs(4326))
  temp = mygeom
  
  extent.nf = extent(bbox(st_coordinates(temp$geometry))) ## better use extent of coastline
  W = owin(c(extent.nf@xmin-1, extent.nf@xmax+1), c(extent.nf@ymin-2, extent.nf@ymax+2))
  centroids.df = st_coordinates(temp$geometry)
  ppp.nf = as.ppp(centroids.df, W=W)

  #' Pairwise correlation function, Ripley
  pcf_result = Kinhom(ppp.nf, correction = "Ripley")
  g <- pcf_result$iso # --> The PCF values with Ripley's isotropic edge correction.
  # Remove Inf and NA values if required
  #valid =  g #<- is.finite(g) & !is.na(g)
  my_list[i,]$Ripley = mean(g, na.rm=T)
}

my_list_tot = my_list %>%left_join(ATot, by="Year")

my_list_tot = my_list_tot%>%left_join(CC_GAMM, by="Year")%>%
  left_join(UA_GAMM, by="Year")%>%
  mutate(Pred_Tot = Pred_UI + (Pred_CC*2))


#' Fit trendline
gam_rip = gam(Ripley ~ s(Year, k=3), family="gaussian", data=my_list_tot)
coeff <- 100
my_list_tot$fit = gam_rip$fitted.values

#' plot Ripley over time (Years)
ggplot(my_list_tot, aes(x=Year, y=Ripley))+
  geom_point(aes(y=Ripley), alpha=1, size=1)+
  geom_line(aes(y=fit), lwd=1, alpha=.8, col="blue")+
  ylab("Ripley's K")+
  theme_bw()

# Values > 1 → overall clustering
# Values < 1 → overall repulsion

