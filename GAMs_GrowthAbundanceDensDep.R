#' R Code accompanying research article
#' Grundlehner et al. (2026) 
#' TITLE
#' JOURNAl
#' DOI:
#'  -- (3)  Fit statistical models for temporal trends of SRW abundances by aggregation area + Test density dependent effects
#' Script by: Anne Grundlehner (anne.grundlehner@utas.edu.au), 
#' Affiliation: Institute of Marine & Antarctic studies, Hobart, Tasmania, Australia
#'
#'#' **** ANY USE OF THIS CODE SHOULD REFER TO GRUNDLEHNER ET AL (2026) ****


#### ------------- INITIATION ---------------- ####

rm(list=ls())

library(readr)
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(mgcv)
library(rnaturalearthdata)
library(rnaturalearth)

#' Set your working directory
my_wd = "~/Documents"
setwd(my_wd) 

setwd("~/Documents/PhD/Ch2/Script/Clean scripts github")


#### ------------ Import data ------------ ####

#' Import Australia outline for plotting
ossie <- ne_countries(scale = "medium", returnclass = "sf", country = "australia")
p <- ggplot(data = ossie) + geom_sf() + xlab("Longitude") + ylab ("Latitude")+
  coord_sf(xlim = c(114, 133), ylim = c(-30, -36))+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,"cm"),
        strip.background = element_rect(fill = "white"), # background col for facet grid
        legend.text = element_text(size = 12, colour = " black"))
p


#' Import abundance data by aggregation area
fg_area = read_csv("fg_area_GH.csv")

#' Import the results from cluster analysis:
my_polys = read_sf("clusteranalysis_polygons_aggregation_areas.shp")%>%
  arrange(-desc(Lon_centre)) #' IMPORTANT


#### ----------- 0. Housekeeping ---------- ####


#' Correct names
my_polys$name = c("Cape leeuwin", "Albany to Pallinup",
                  "Bremer and Doubtful Islands", "Cape le Grand",
                  "Cape arid", "Recherche and Israelite Bay",
                  "Twilight cove", "Eucla", "Nullarbor", "Fowlers Bay")
poly_labels = c("Cape Leeuwin", "Albany to Cheynes Bay",
                  "Bremer to Doubtful Islands Bays", "Cape le Grand",
                  "Cape Arid", "Israelite Bay to Pt Culver",
                  "Twilight Cove", "Eucla", "Nullarbor", "Fowlers Bay")
my_polys$name = factor(my_polys$name, levels=c("Cape leeuwin", "Albany to Pallinup",
                                               "Bremer and Doubtful Islands", "Cape le Grand",
                                               "Cape arid", "Recherche and Israelite Bay",
                                               "Twilight cove", "Eucla", "Nullarbor", "Fowlers Bay"),
                       labels = poly_labels)

# Add aggregation area centre
my_poly_info = my_polys%>%
  as.data.frame()%>%dplyr::select(name, Lon_centre, Lat_centre)%>%
  rename(aggregation_area = "name")

fg_area = fg_area %>% left_join(my_poly_info, by="aggregation_area")

area_arranged = my_poly_info %>%
  dplyr::select(aggregation_area, Lon_centre)%>%
  arrange(-desc(Lon_centre))%>%
  dplyr::select(aggregation_area)%>%
  mutate(aggregation_area = as.character(aggregation_area))%>%
  rbind(c("OTHER"))

area_arranged$aggregation_area = factor(area_arranged$aggregation_area, levels=c("Cape leeuwin", "Albany to Pallinup",
                                               "Bremer and Doubtful Islands", "Cape le Grand",
                                               "Cape arid", "Recherche and Israelite Bay",
                                               "Twilight cove", "Eucla", "Nullarbor", "Fowlers Bay", "OTHER"),
                       labels = c("Cape leeuwin", "Albany to Pallinup",
                                  "Bremer and Doubtful Islands", "Cape le Grand",
                                  "Cape arid", "Recherche and Israelite Bay",
                                  "Twilight cove", "Eucla", "Nullarbor", "Fowlers Bay", "OTHER"))


#### ------- 1. Retrieve lag data -------  ####

#' Loop through areas to retrieve 'Initial' Abundances 

#' Retrieve 'Initial' Abundances 
#' 5 = first 5 years of the study (1977-1982)
#' 10 = first 10 years of the study  (1977-1986)
#' 20 = first 20 years of the study (1977-1996)
fg_area$CalfInit5 = NA
fg_area$TotInit5 = NA
fg_area$UnaccInit5 = NA
fg_area$CalfInit10 = NA
fg_area$TotInit10 = NA
fg_area$UnaccInit10 = NA
fg_area$CalfInit20 = NA
fg_area$TotInit20 = NA
fg_area$UnaccInit20 = NA


for(i in 1:nrow(fg_area)){
  y = fg_area[i,]$Year
  loc = fg_area[i,]$aggregation_area
  dat_init20 = fg_area[fg_area$aggregation_area == loc & fg_area$Year %in%c(1977:1996)==T,]
  dat_init10 = fg_area[fg_area$aggregation_area == loc & fg_area$Year %in%c(1977:1986)==T,]
  dat_init5 = fg_area[fg_area$aggregation_area == loc & fg_area$Year %in%c(1977:1982)==T,]
  dat = fg_area[fg_area$aggregation_area == loc & fg_area$Year %in%c(y, y-1, y-2)==T,]
  dat2 = fg_area[fg_area$aggregation_area == loc & fg_area$Year %in%c(y-1, y-2, y-3)==T,]

  #' Add initial occuppance as a variable as well
  if(nrow(dat_init5)>=1){
    fg_area[i,]$CalfInit5 = mean(dat_init5$calf, na.rm=T)
    fg_area[i,]$TotInit5 = mean(dat_init5$total, na.rm=T)
    fg_area[i,]$UnaccInit5 = mean(dat_init5$unacc, na.rm=T)
  }else{
    fg_area[i,]$CalfInit5 = NA
    fg_area[i,]$TotInit5 = NA
    fg_area[i,]$UnaccInit5 = NA
  }
  if(nrow(dat_init10)>=1){
    fg_area[i,]$CalfInit10 = mean(dat_init10$calf, na.rm=T)
    fg_area[i,]$TotInit10 = mean(dat_init10$total, na.rm=T)
    fg_area[i,]$UnaccInit10 = mean(dat_init10$unacc, na.rm=T)
  }else{
    fg_area[i,]$CalfInit10 = NA
    fg_area[i,]$TotInit10 = NA
    fg_area[i,]$UnaccInit10 = NA
  }
  if(nrow(dat_init20)>=1){
    fg_area[i,]$CalfInit20 = max(dat_init20$calf, na.rm=T)
    fg_area[i,]$TotInit20 = max(dat_init20$total, na.rm=T)
    fg_area[i,]$UnaccInit20 = max(dat_init20$unacc, na.rm=T)
  }else{
    fg_area[i,]$CalfInit20 = NA
    fg_area[i,]$TotInit20 = NA
    fg_area[i,]$UnaccInit20 = NA
  }
  
}


#' Add lon lat back to area
fg_area = fg_area %>% left_join(my_poly_info, by="aggregation_area")




#### ------------------- 2.1 BASE MODEL: ABUNDANCE -----------------------####

fg_area = fg_area%>%dplyr::filter(aggregation_area != "OTHER")

g0 = gam(calf ~ 1 + s(Year),
         family="poisson", data = fg_area)

fg_area$residuals = g0$residuals

#' check for patterns in residuals
ggplot(fg_area, aes(x=calf, y=residuals))+
  geom_point()+
  facet_wrap(.~aggregation_area, scales="free")+
  geom_smooth(method="gam", formula=y~s(x, k=3))
# Conclusion: remaining patterns by area

#' add area as smoother interaction
g0 = gam(calf ~  s(Year, k=3, by=factor(aggregation_area)) + factor(aggregation_area),
         family="poisson", data = fg_area)



#### ------------------- 2.2. MODEL VALIDATION -----------------------####

#' check for patterns in residuals
par(mfrow=c(2,2))
gam.check(g0)
par(mfrow=c(1,1))


#' Plot this model
pred_data <- expand_grid(
  Year = seq(from=min(fg_area$Year), 
             to=max(fg_area$Year), 
             length.out = 10000),
  aggregation_area = levels(fg_area$aggregation_area))
dim(pred_data)

pred_data <- pred_data %>% 
  predict(g0, newdata = ., se=TRUE) %>% 
  as_tibble() %>% 
  cbind(pred_data)%>%
  mutate(SeUp = fit + se.fit,
         SeLo = fit - se.fit)

#' Model validation G0
plot(g0$fitted.values ~ fg_area$calf)
hist(g0$fitted.values) #  GAM --> Homoscedasticity

ggplot(fg_area)+
  geom_point(aes(y=g0$residuals, x=g0$fitted.values))+
  facet_wrap(~aggregation_area, scales="free")+
  geom_hline(yintercept=0, col="red") # residuals vs fitted

ggplot(fg_area)+
  geom_point(aes(y=g0$residuals, x=calf))+
  facet_wrap(~aggregation_area, scales="free") # residuals vs observed

ggplot(fg_area)+
  geom_point(aes(y=g0$residuals, x=Year))+
  facet_wrap(~aggregation_area, scales="free")+
  geom_hline(yintercept=0, col="red")+ # residuals ~ year
  geom_smooth(aes(y=g0$residuals, x=Year))
#'acceptable

#' Plot the fitted value for the temporal trend (CALF)
fg_area$F1 = g0$fitted.values
p <- ggplot()
p <- p + geom_point(data = fg_area, 
                    aes(y = F1, x = Year),
                    shape = 16, 
                    size = 2)
p <- p + geom_hline(yintercept = 0) + ylab("Fitted values")+xlab("Years")
p <- p + theme(text = element_text(size=15))
p <- p + facet_wrap( ~ aggregation_area, scales = "fixed")
p # smoothers per route

#' Plot residuals
fg_area$E1 = g0$residuals
p <- ggplot()
p <- p + geom_point(data = fg_area, 
                    aes(y = E1, x = Year),
                    shape = 16, 
                    size = 2)
p <- p + geom_hline(yintercept = 0) + ylab("Residuals")+xlab("Years")
p <- p + theme(text = element_text(size=15))
p <- p + facet_wrap( ~ aggregation_area, scales = "fixed")
p # Residuals per route: ok. maybe remove outliers Aba-Palli????


# smoothers per route
p <- ggplot(pred_data, aes(col=aggregation_area))+
  #geom_point(data=fg_area, aes(y=calf, x=Year), alpha=0.4)+
  geom_line(aes(y = fit, x = Year), size = 1)+
  geom_line(aes(y = SeLo, x = Year), size = 0.5, lty="dashed")+
  geom_line(aes(y = SeUp, x = Year), size = 0.5, lty="dashed")+
  ylab("Fitted values")+xlab("Year")+
  theme(text = element_text(size=15))+ 
  facet_wrap( ~ aggregation_area, scales = "free_y")
p 



#### ------------------- 2.3 Plot final BASE MODEL -----------------------####

fg_area$F1 = g0$fitted.values

p <- ggplot(fg_area, aes(col=aggregation_area))+
  geom_line(aes(y = F1, x = Year), size = 1)+
  geom_point(aes(y=calf, x=Year), alpha=0.4)+
  ylab("Abundance of CC pairs (n pairs)")+xlab("Year")+
  theme(text = element_text(size=15), legend.position = "none")+ 
  facet_wrap( ~ aggregation_area, scales = "free_y", ncol=2)+
  theme_bw()+labs(col="Area")
p # smoothers per route


#### ---> PLOT TEMPTREND  MODELS ALL DEMO GROUPS IN 1 ####
temp1 = fg_area %>% dplyr::select(Year, aggregation_area, calf)%>%
  rename(Count = calf)%>%mutate(Class = "CC pairs")%>%
  mutate(aggregation_area = factor(aggregation_area))
temp2 = fg_area %>% dplyr::select(Year, aggregation_area, unacc)%>%
  rename(Count = unacc)%>%mutate(Class = "UA")%>%
  mutate(aggregation_area = factor(aggregation_area))
temp3 = fg_area %>% dplyr::select(Year, aggregation_area, total)%>%
  rename(Count = total)%>%mutate(Class = "Total")%>%
  mutate(aggregation_area = factor(aggregation_area))

g1 = gam(Count ~ 1 +s(Year, by = aggregation_area, k=3)+aggregation_area,
         family="poisson", data = temp1)
g2 = gam(Count ~ 1 + s(Year, by=aggregation_area,  k=3)+aggregation_area,
         family="poisson",data = temp2)
g3 = gam(Count ~ 1 + s(Year, by=aggregation_area, k=3)+aggregation_area,
         family="poisson",data = temp3)

temp1$F1 = g1$fitted.values
temp2$F1 = g2$fitted.values
temp3$F1 = g3$fitted.values
temp = rbind(temp1, temp2, temp3)
temp$Class = factor(temp$Class, levels = c("CC pairs", "UA", "Total"),
                    labels = c("Cow-calf pairs", "Unaccompanied individuals", "All individuals"))

p <- ggplot(temp, aes(col=Class))+
  geom_point(aes(y=Count, x=Year), alpha=0.3, size=0.8)+
  geom_line(aes(y = F1, x = Year), size = 1)+
  ylab("Number of individuals")+xlab("Year")+
  scale_color_manual(values = c("coral2", "turquoise3", "grey70"))+
  facet_wrap( ~ aggregation_area, scales = "free_y", ncol=2)+
  theme_bw()+theme(text = element_text(size=10), legend.position = "bottom")+ labs(col="Group")
p 



MC = c(brewer.pal(n=8, "Dark2")[1:6], brewer.pal(n=12, "Paired")[c(2,5,12,8)])
p <- ggplot(temp, aes(col=aggregation_area))+
  #geom_point(aes(y=Count, x=Year), alpha=0.3, size=0.8)+
  geom_line(aes(y = F1, x = Year), size = 1)+
  ylab("Abundance (n)")+xlab("Year")+
  scale_color_manual(values = MC)+
  facet_wrap(~Class, scales = "free_y", ncol=1)+
  theme_bw()+labs(col="Aggregation area")+
  theme(text = element_text(size=13), 
        legend.title = element_text(face="bold"))
  
p1=p



#### ---------- 3.  Plot Growth based on fitted values ------------- ####

temp1 = fg_area %>% dplyr::select(Year, aggregation_area, calf)%>%rename(Count = calf)%>%mutate(Class = "CC pairs")
temp2 = fg_area %>% dplyr::select(Year, aggregation_area, unacc)%>%rename(Count = unacc)%>%mutate(Class = "UA")
temp3 = fg_area %>% dplyr::select(Year, aggregation_area, total)%>%rename(Count = total)%>%mutate(Class = "Total")

g1 = gam(Count ~ 1 + s(Year, by=as.factor(aggregation_area), k=3) + aggregation_area,
         family="poisson", data = temp1)
g2 = gam(Count ~ 1 + s(Year, by=as.factor(aggregation_area),  k=3)+aggregation_area,
         family="poisson",data = temp2)
g3 = gam(Count ~ 1 + s(Year, by=as.factor(aggregation_area), k=3)+aggregation_area,
         family="poisson",data = temp3)

temp1$F1 = g1$fitted.values
temp2$F1 = g2$fitted.values
temp3$F1 = g3$fitted.values
temp = rbind(temp1, temp2, temp3)
temp$Class = factor(temp$Class, levels = c("CC pairs", "UA", "Total"),
                    labels = c("Cow-calf pairs", "Unaccompanied individuals", "All individuals"))

temp$Growth = NA 
for(i in 1:nrow(temp)){
  j = temp[i,]$Year
  a = temp[i,]$aggregation_area
  c= temp[i,]$Class
  jm1 = j-1
  sub = temp[temp$Year==jm1 & temp$aggregation_area==a & temp$Class == c,]
  if(nrow(sub)>0){
    temp[i,]$Growth = (temp[i,]$F1 - sub$F1) / sub$F1
  }
}
summary(temp$Growth)


#### --> PLOT GROWTH AND ABUNDANCE FOR FIGURE USED IN MANUSCRPT ####
p2 = temp %>%
  ggplot(aes(x=Year, y=Growth*100, 
             col=aggregation_area, group=aggregation_area))+
  #geom_point()+
  geom_smooth()+
  facet_grid(Class~.)+
  ylab("Relative rate of increase (%)")+
  geom_hline(yintercept = 0, lty="dashed")+
  scale_color_manual(values=MC)+
  ylim(-10,50)+
  #theme(text = element_text(size=15),
  #      axis.text  = element_text(size=15),)+
  theme_bw()
p2

#' ABUNDANCE
p <- ggplot(temp, aes(col=aggregation_area))+
  theme_bw()+
  geom_hline(yintercept = 0, col="grey", lty="dotted")+
  #geom_point(aes(y=CAGR, x=Year), alpha=0.3, size=1)+
  scale_colour_manual(values=MC)+
  geom_line(aes(y = F1, x = Year), size = 0.7)+
  ylab("Annual abundance (n)")+xlab("Year")+
  theme(#text = element_text(size=15), 
        #axis.text  = element_text(size=15),
    legend.key.spacing.y = unit(0, "pt"),, 
    legend.key.spacing.x = unit(2, "pt"),,
        legend.title = element_text(face="bold"))+
  guides(colour = guide_legend(nrow = 4)) +
  facet_grid(Class~., scales = "free_y")+
  labs(col="Aggregation area")
p
p1=p

ggarrange(p1, p2, common.legend = T, 
          legend="bottom")






#### --------------4. Investigate effect of Initial Abundance -----------------------####

#' Use base model --> but exclude observations from the first 10 or 20 years (bias otherwise)

g0 = gam(calf ~ 1+  s(Year, by=factor(aggregation_area), k=3) + aggregation_area,
         family="poisson", data = fg_area[fg_area$YearBin != "1976-1985",])

g1b = gam(calf ~ 1+ CalfInit10 + s(Year, by=factor(aggregation_area), k=3)+ aggregation_area,
          family="poisson",data = fg_area[fg_area$YearBin != "1976-1985",])
summary(g1b)
AIC(g0, g1b)

g0 = gam(calf ~ 1+  s(Year, by=factor(aggregation_area), k=3),
         family="poisson",data = fg_area[fg_area2$YearBin %in%c("1976-1985", "1986-1995") ==F,])
summary(g0)

g1c = gam(calf ~ 1+ CalfInit20 + s(Year, by=factor(aggregation_area), k=3) + 
            aggregation_area,
          family="poisson",data = fg_area[fg_area$YearBin %in%c("1976-1985", "1986-1995") ==F,])
summary(g1c)

AIC(g0, g1c)

#' does the effect hold up for all periods?
g1b = gam(calf ~ 1+ CalfInit10*YearBin + s(Year, by=factor(aggregation_area), k=3)+ aggregation_area,
family="poisson",data = fg_area[fg_area$YearBin != "1976-1985",])
summary(g1b)

g1c = gam(calf ~ 1+ CalfInit20*YearBin + s(Year, by=factor(aggregation_area), k=3)+ 
            aggregation_area,
          family="poisson",data = fg_area[fg_area$YearBin %in%c("1976-1985", "1986-1995") ==F,])
summary(g1c)

#' total
gt = gam(total ~ 1+ TotInit10*YearBin + s(Year, by=factor(aggregation_area), k=3)+ 
            aggregation_area,
          family="poisson",data = fg_area[fg_area$YearBin != "1976-1985",])
summary(gt)

#' ua
g1c = gam(total ~ 1+ UnaccInit10*YearBin + s(Year, by=factor(aggregation_area), k=3)+ 
 aggregation_area,family="poisson",data = fg_area[fg_area$YearBin != "1976-1985",])
summary(g1c)





#### ------ 5. Density depedent effects: Do  local abundance respond to regional abundances? ------- ####


# Calf ~ calf
g0 = g1 = gam(calf ~ 1 + s(Year, k=4, by=factor(aggregation_area)), 
              family="poisson",data = fg_area)
g1 = gam(calf ~ 1 + CalfYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)

# Calf ~ tot
g1 = gam(calf ~ 1 + TotYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)

# Unacc ~ unacc
g0 = g1 = gam(unacc ~ 1 + s(Year, k=4, by=factor(aggregation_area)), family="poisson",data = fg_area)
g1 = gam(unacc ~ 1 + UnaccYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)

# Unacc ~ tot
g1 = gam(unacc ~ 1 + TotalYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)


# Tot ~ calf
g0 = g1 = gam(total ~ 1 + s(Year, k=4, by=factor(aggregation_area)), family="poisson",data = fg_area)
g1 = gam(total ~ 1 + CalfYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)

# Tot ~ unacc
g0 = g1 = gam(total ~ 1 + s(Year, k=4, by=factor(aggregation_area)), family="poisson",data = fg_area)
g1 = gam(total ~ 1 + UnaccYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)

# Tot ~ tot
g1 = gam(total ~ 1 + TotalYearSum*aggregation_area +
           s(Year, k=4, by=factor(aggregation_area)) ,
         family="poisson", data = fg_area)
summary(g1)
AIC(g0, g1)



#### ------ 6. Density depedent effects: Do local abundances respond to local abundances of other demographic groups? ------- ####

fg_area2 = fg_area2

#### --------------- 6.1. Local spatial segregation?  -------------------------- ####

ggplot(fg_area2, aes(x=calf, y=unacc))+geom_point()+
  geom_smooth(method="gam", formula=y~s(x,k=3))+xlab("Local CC abundance")+
  facet_wrap(~aggregation_area, scales="free", ncol=5)

#' is there a negative relation between CC and UA
#' account for the overall year effect
g1 = gam(unacc ~ 1 + aggregation_area*calf, family="poisson", data=fg_area)
summary(g1)

#' OR: is the variation from the general trend explained by CC ?
g0 = gam(unacc ~ 1 + s(Year, by=factor(aggregation_area)) , family="poisson", data=fg_area)
fg_area$FU = g0$fitted.values
#' Does variation from the mean trend negatively relate to CC ?
g1 = gam(unacc ~ 1 + FU + calf, family="poisson", data=fg_area)
summary(g1)

#' No, there is no segregation / UA locally driven away.



#### --------------- 6.2 Selective use of areas depending on total annual SRW abundance?  -------------------------- ####


ggplot(fg_area, aes(y=total, x=TotalYearSum))+geom_point()+
  geom_smooth(method="gam", formula=y~s(x,k=3))+xlab("Local abundance")+
  facet_wrap(~aggregation_area, scales="fixed", ncol=5)

#' Do certain areas respond differently to TotYearSum than others? (A)
g0 = gam(total ~ 1 + s(Year, by=factor(aggregation_area)) + 
           TotalYearSum + aggregation_area, data=fg_area, 
         family="gaussian")
g1 = gam(total ~ 1 + s(Year, by=factor(aggregation_area)) + 
           TotalYearSum * aggregation_area, data=fg_area, 
         family="gaussian")
AIC(g0, g1)
summary(g1)

#' And if there is locally a massive peak or dip, is that explained variation somewhere? (B)
#' Aka is there a relation between the residuals over time
g0 = gam(total ~ 1 + s(Year, by=factor(aggregation_area), k=3), data=fg_area, 
         family="poisson")
fg_area$Res =  g0$fitted.values - fg_area$total



#' Repeat for CC

#' Do certain areas respond differently to TotYearSum than others? (A)
g0 = gam(calf ~ 1 + s(Year, by=factor(aggregation_area), k=5) + 
           CalfYearSum + aggregation_area, data=fg_area, 
         family="gaussian")
g1 = gam(calf ~ 1 + s(Year, by=factor(aggregation_area), k=5) + 
           CalfYearSum * aggregation_area, data=fg_area, 
         family="gaussian")
AIC(g0, g1)
summary(g1)




#### ---------- 6.2 response of calf to  calfyearsum (not annual regional total, but regional CC total---------------  ####
 
for(i in 1:length(unique(fg_area$aggregation_area))){
  print("------------------------------")
  print(unique(fg_area$aggregation_area)[i])
  
  g0 = gam(calf ~ 1 + s(Year) + CalfYearSum, 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], 
           family="poisson")
  print(summary(g0))
  g1 = gam(calf ~ 1 + s(Year)+ s(CalfYearSum, k=5), 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], family="poisson")
  if(AIC(g0) > (AIC(g1)+10)){
    print("** non-linear **")
    print(AIC(g0, g1))
  }
}


#### ---------- 6.3 response of CC to YearSum (regional total) ---------------  ####
for(i in 1:length(unique(fg_area$aggregation_area))){
  print(unique(fg_area$aggregation_area)[i])
  g0 = gam(calf ~ 1 + s(Year) + TotalYearSum, 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], 
           family="poisson")
  print(summary(g0))
  g1 = gam(calf ~ 1 + s(Year)+ s(TotalYearSum, k=5), 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], family="poisson")
  if(AIC(g0) > (AIC(g1)+10)){
    print("** non-linear **")
    print(AIC(g0, g1))
  }
}

#### ---------- 6.4 response of UA to totalyearsum (regional total) ---------------  ####
#' Is there evidence of UA being driven away from certain areas (i.e. as proposed by O'shannessey et al 2025)
for(i in 1:length(unique(fg_area$aggregation_area))){
  print("------------------------------")
  print(unique(fg_area$aggregation_area)[i])
  g0 = gam(unacc ~ 1 + s(Year) + TotalYearSum, 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], family="poisson")
  print(summary(g0))
  g1 = gam(unacc ~ 1 + s(Year)+ s(TotalYearSum, k=5), 
           data=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ], family="poisson")
  if(AIC(g0) > (AIC(g1)+10)){
    print("** non-linear **")
    print(AIC(g0, g1))
  }
}


#' Again, but force fitted values to get pure calf effect.
for(i in 1:length(unique(fg_area$aggregation_area))){
  print("------------------------------")
  print(unique(fg_area$aggregation_area)[i])
  sub=fg_area[fg_area$aggregation_area ==unique(fg_area$aggregation_area)[i],  ]
  gx = gam(unacc ~ 1 + s(Year), 
           data=sub, family="poisson")
  sub$F1 = gx$fitted.values
  g0 = gam(unacc ~ 1 + F1 + TotalYearSum, 
           data=sub, family="poisson")
  if(AIC(gx) > (AIC(g0)+10)){
    print("** important **")
    print(AIC(gx, g0))
  }
  print(summary(g0))
  g1 = gam(unacc ~ 1 + F1 + s(TotalYearSum, k=5), 
           data=sub, family="poisson")
  if(AIC(g0) > (AIC(g1)+10)){
    print("** non-linear **")
    print(AIC(g0, g1))
  }
}










