#' R Code accompanying research article
#' Grundlehner et al. (2026) 
#' TITLE
#' JOURNAl
#' DOI:
#'  -- (2)  Create aggregation areas (2-step Cluster analysis)
#' Script by: Anne Grundlehner (anne.grundlehner@utas.edu.au), 
#' Affiliation: Institute of Marine & Antarctic studies, Hobart, Tasmania, Australia
#' 
#' **** ANY USE OF THIS CODE SHOULD REFER TO GRUNDLEHNER ET AL (2026) ****


rm(list=ls())

#' Set your working directory
my_wd = "~/Documents"
setwd(my_wd) 

#' load libraries
library(sp)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(igraph)
library(alphahull)
library(rgeos)
library(spatstat)
library(rnaturalearth)
library(rnaturalearthdata)
library(xlsx)
library(readlxl)
library(readxl)


#' Simulate distance between cornering gridcels
p1 = c(120, -35)
p2 = c(120.1, -35.1)
p1 = st_sfc(st_point(p1))%>% st_set_crs(4326)
p2 = st_sfc(st_point(p2))%>% st_set_crs(4326)
st_distance(p1, p2)/1000
dist_to_cornering = st_distance(p1, p2)/1000
dist_to_cornering 

#' Simulate distance between neighboring grid cells
p1 = c(120.1, -35)
p2 = c(120, -35)
p1 = st_sfc(st_point(p1))%>% st_set_crs(4326)
p2 = st_sfc(st_point(p2))%>% st_set_crs(4326)
st_distance(p1, p2)/1000
dist_to_neighbouring_x = st_distance(p1, p2)/1000
dist_to_neighbouring_x

p1 = c(120, -35)
p2 = c(120, -35.1)
p1 = st_sfc(st_point(p1))%>% st_set_crs(4326)
p2 = st_sfc(st_point(p2))%>% st_set_crs(4326)
st_distance(p1, p2)/1000
dist_to_neighbouring_y = st_distance(p1, p2)/1000
dist_to_neighbouring_y

#' Note the difference between x and y! 
#' Most of our grid cells are aligned on the y-axis
#' To include neighboring grid cells, set dist to < 12
#' To include cornering grid cells, set dist to < 15



#### INITIATION: prepare functions ####
# function to convert an arc into line segments
# given the center of the arc, the radius, the vector, and the angle (radians)
arc2line <- function(center, r, vector, theta, npoints = 100) {
  angles <- anglesArc(vector, theta)
  seqang <- seq(angles[1], angles[2], length = npoints)
  x <- center[1] + r * cos(seqang)
  y <- center[2] + r * sin(seqang)
  coords.xy <- cbind(x,y)
  line <- Line(coords = coords.xy)
  return(line)
}


#' Hulls to lines:
ahull2lines <- function(hull){
  arclist <- hull$arcs
  lines <- list()
  for (i in 1:nrow(arclist)) {
    # Extract the attributes of arc i
    center_i <- arclist[i, 1:2]
    radius_i <- arclist[i, 3]
    vector_i <- arclist[i, 4:5]
    theta_i <- arclist[i, 6]
    # Convert arc i into a Line object
    line_i <- arc2line(center = center_i, r = radius_i, vector = vector_i, theta = theta_i)
    list_length <- length(lines)
    if(list_length > 0){
      last_line_coords <- lines[[list_length]]@coords
    }
    if(i == 1){
      # Add the first line to the list of lines
      lines[[i]] <- line_i
    } else if(isTRUE(all.equal(line_i@coords[1,], last_line_coords[nrow(last_line_coords),]))){
      lines[[list_length]]@coords <- rbind(last_line_coords, line_i@coords[2:nrow(line_i@coords),])
    } else {
      lines[[length(lines) + 1]] <- line_i
    }
  }
  # Convert the list of lines to a Line object
  lines <- Lines(lines, ID = 'l')
  # Convert the Line object to a SpatialLines object
  sp_lines <- SpatialLines(list(lines))
  return(sp_lines)
}


#' convert the SpatialLines object into a SpatialPolygon
spLines2poly <- function(sp_lines){
  lines_slot <- sp_lines@lines[[1]]
  poly_bool <- sapply(lines_slot@Lines, function(x){
    coords <- lines_slot@Lines[[1]]@coords
    all.equal(coords[1,], coords[nrow(coords),])
  })
  poly_lines <- sp_lines[poly_bool]
  poly_lines_slot <- poly_lines@lines
  sp_polys <- SpatialPolygons(list(Polygons(lapply(poly_lines_slot, function(x) {
    Polygon(slot(slot(x, "Lines")[[1]], "coords"))
  }), ID = "1")))
  return(sp_polys)
}


#' Combined: alpha hull to spatial polygons
ahull2poly <- function(hull){
  # Convert the alpha hull to SpatialLines
  hull2SpatialLines <- ahull2lines(hull)
  # Convert SpatialLines to SpatialPolygon
  SpatialLines2SpatialPolygon <- spLines2poly(hull2SpatialLines)
  return(SpatialLines2SpatialPolygon)
}



####  INITIATION - Prepare data and plots bases ####

#' Import australia
ossie <- ne_countries(scale = "medium", returnclass = "sf", country = "australia")
p <- ggplot(data = ossie) + geom_sf() + xlab("Longitude") + ylab ("Latitude")+
  coord_sf(xlim = c(114, 133), ylim = c(-30, -36))+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,"cm"),
        strip.background = element_rect(fill = "white"), # background col for facet grid
        legend.text = element_text(size = 12, colour = " black"))
p




#' Import Whale data

# Gridded data, summarized per gridcel per year (i.e. no duplicates, but includes zero obs.)
df = read_excel("fg_grid0.1_AggrAreaAttributed.xlsx")%>% 
  mutate(CC = calf*2, # because CC pair = 2 indiv
         unacc = ifelse(unacc < 0, 0, unacc))%>%
  dplyr::select(!c(calf, total))%>%
  pivot_longer(cols = c(unacc, CC), 
               names_to = "Class", values_to="Count")%>%
  mutate(Lon_decimal = LonBin,
         Lat_decimal = LatBin)%>%
  #' we are not interested in zero observations here, so drop that:
  dplyr::filter(Count > 0)

summary(df$Count) #' Check
temp = df[df$Count > 100,]

df_cluster = df %>%
  arrange(-desc(Year))

#' make sf object
mygeom=(st_as_sf(x=df_cluster, 
                 coords = c("Lon_decimal", "Lat_decimal"), 
                 crs=crs(ossie)))

#### ----------- 0. Set parameter preferences --------- ####


#' Maximum distance between observations within 1 cluster
max_dist = 15 #' (note: distance in km)

#' Do you want to plot the clusters whilst looping?
plot_cluster = T

#' make minimum number of whales in cluster dependent on year
minsrw_dep = T

#' do you want to base the cluster analysis on CC only or also on UA
calf_only = F # calf only

# because non-cluster data is removed; and the ACTUAL occupation is checked later.
minimum_years = 5 # The (combined) cluster must have 

# AND/OR use by X animals? (not used)
minimum_cc_alltime = 5

cluster_period = "NA"

#' Set Maxdist for the combination for the second loop, where we combine years
#' we are only looking for overlap BETWEEN years
#' Just let clusters that are neighbouring and cornering be combined
#' but if there is more than [max_dist2] km in between, AND they are in different years,
#' then it does not make sense to combine them;
#' Note that we overlay ALL years; so if this threshold is not reached,
#'  they are quite consistently seperated. 
max_dist2 = 15


#### ------- Part 1: PLOTTING + MAKE CLUSTERS BY YEAR -- Vary max_dist (yearly clusters) --------------------- ####


#### ------------------ 1.1 Cluster analyses (part 1) --------------------- ####

all_clust_res = data.frame()

#' Loop through cluster periods and run cluster analysis to identify
#' locations where whale aggregations are found, for each year.
for(i in unique(df_cluster$Year)){

  mydat0 = mygeom[mygeom$Year == i,] # Year
  mydat0$period_id = paste("Y", mydat0$Year, sep="") #' Set period ID

  #' now filter the class you want to use
  mydat0 = mydat0 %>% 
    mutate(Class = ifelse(Class == "unacc", "UA", Class))
  if(calf_only == T){
    mydat = mydat0%>%dplyr::filter(Class == "CC")
  }else{
    mydat = mydat0
  }
  
  #' Check:
  unique(mydat$Class)
  nrow(mydat)
  
  if(nrow(mydat) > 1){ # If there is only 1 obs in the entire year, don't run cluster analysis.
    
    distmat = st_distance(mydat$geometry) 
    distmat = distmat/1000 # covert to km
    
    #' run cluster analysis on spatial points where observations were made
    chc <- hclust(as.dist(distmat), method="single")
    chc.d40 <- cutree(chc, h=max_dist) 
    xy = as(mydat$geometry, "Spatial")
    
    # Join results 
    xy$data <- data.frame(xy, Clust=chc.d40)
    xy@data <- data.frame(xy, Clust=chc.d40)
    xy = as.data.frame(xy)
    
    #' Bind xy back to count data
    xy = xy%>%
      dplyr::rename(LonBin = "coords.x1", LatBin = "coords.x2")
    nrow(xy)
    
    xy <- xy %>% dplyr::select(LonBin, LatBin, Clust)%>%
      left_join(mydat, by=c("LonBin", "LatBin"))
    xy =   as.data.frame(lapply(xy, rep, xy$Count)) # longlist for each obs = 1 point
    
    #' Remove clusters with < X whales
    #' If we made it dependent on the total number of whales, 
    #' first calculate min_srw:
    if(minsrw_dep == T){
      if(unique(xy$Year) < 1990){
        min_srw_cluster = 1
      }else{
        if(unique(xy$Year) < 2001){
          min_srw_cluster = 2
        }else{
          if(unique(xy$Year) < 2010){
            min_srw_cluster = 3
          }else{min_srw_cluster = 3}
        }}}
    
    #' MAKE SURE IT IS A LONGLIST with each point being a single whale
    for(i in unique(xy$Clust)){
      n = nrow(xy[xy$Clust==i,])
      if(n < min_srw_cluster){
        xy[xy$Clust==i,]$Clust = "X"
      }}
    
    #plot_cluster = T
    if(plot_cluster == T){
      myplot = ggplot(xy[,])+
        geom_sf(data=ossie$geometry)+
        geom_point(alpha=0.5,
                   aes(x=LonBin, y=LatBin, col=factor(Clust)))+
        #scale_colour_viridis_d(option="F")+
        coord_sf(xlim = c(114, 133), ylim = c(-30, -36))+
        ggtitle(xy$Year)
      print(myplot)
    }
    
    my_cluster_result = xy
    
    #' Add ID info
    my_cluster_result$cluster_period = cluster_period
    my_cluster_result$max_dist_used = max_dist
    my_cluster_result$min_srw_in_cluster_used = min_srw_cluster
    
    #' Export / Save this round
    if(nrow(all_clust_res) < 1){
      all_clust_res = my_cluster_result
    }else{
      all_clust_res = rbind(all_clust_res, my_cluster_result)
    }
  } 
  } # end loop



#' Check output
head(all_clust_res)
nrow(all_clust_res)
unique(all_clust_res$Year)
unique(all_clust_res$ClustID)
unique(all_clust_res[all_clust_res$Clust != "X",]$Year)


#' retrieve number of clusters for each year:
all_clust_res[all_clust_res$Clust != "X",]%>%
  dplyr::select(Clust, Year)%>%
  distinct()%>%
  mutate(n = 1)%>%
  group_by(Year)%>%dplyr::summarise(n=sum(n))%>% print(n=50)

#' Now, we have FILTERED clusters for each year - but just the grid cells included


####--------------- 1.2. Collide all years -------------------- ####

#' Now, we want to collide all overlapping clusters from all the years
# Note that xy still contains all the POINTS and the cluster label
# and the geometry is the point geometry
# so we just use the cluster label to group and filter out scarce ones
# Filter out all X clusters; do not meet criteria

all_clust_res_firstround = all_clust_res

mydat = all_clust_res%>%
  dplyr::filter(Clust != "X")%>%
  mutate(Clust1 = paste(Year, Clust , sep="_C"))%>%
  dplyr::select(!Clust)%>%
  distinct()

#' Run cluster analysis
distmat = st_distance(mydat$geometry) # geometry is still just the gridcel centre
distmat = distmat/1000 # covert to km -- Takes a while

chc <- hclust(as.dist(distmat), method="single")
chc.d40 <- cutree(chc, h=max_dist2) 
xy = as(mydat$geometry, "Spatial")

# Join results 
xy$data <- data.frame(xy, Clust=chc.d40)
#xy@data <- data.frame(xy, Clust=chc.d40)
xy = as.data.frame(xy)%>%
  rename(LonBin = "coords.x1", 
         LatBin = "coords.x2", 
         Clust = "data.Clust")

#' Bind xy back to count data
xy <- xy %>% dplyr::select(LonBin, LatBin, Clust )%>%
  full_join(mydat, by=c("LonBin", "LatBin"), relationship= "many-to-many")

#' Remove clusters with < X whales
xy$nyears = 0

for(i in unique(xy$Clust)){
  n = length(unique(xy[xy$Clust==i,]$Year))
  xy[xy$Clust==i,]$nyears = n
  if(n < minimum_years){
    xy[xy$Clust==i,]$Clust = "X"
  }}

#' filter only unique points for plotting
xy2 = xy[xy$Clust != "X",] %>%dplyr::select("LonBin", "LatBin", "Clust", "nyears")%>%distinct()

myplot = ggplot(xy2)+
  geom_sf(data=ossie$geometry)+theme_bw()+
  geom_point(alpha=0.8,
             aes(x=LonBin, y=LatBin, col=factor(Clust)))+
  #scale_colour_viridis_d(option="F")+
  coord_sf(xlim = c(114, 133), ylim = c(-30, -36))+
  ggtitle(paste("dist =" , paste(max_dist2, "km", sep=" ")))+
  labs(col="Cluster ID")
myplot

xy2$LonBinX = xy2$LonBin
xy2$LatBinX = xy2$LatBin
xy3 <- st_as_sf(xy2, coords=c("LonBinX", "LatBinX"), crs=crs(x))

for(i in unique(xy2$Clust)){
  print(i)
  hulldat =  xy2[xy2$Clust == i,]
  hulldat = hulldat %>%dplyr::arrange(desc(LonBin))%>% arrange(LatBin)%>%
    dplyr::select(LonBin, LatBin)%>%distinct()
  
  #' note that, if the structure is linear, ahull does not work.
  #' Work around that by creating a buffer zone manually instead of using ahull
  #' 1. Set up a test:
  if(length(unique(hulldat$LonBin))==1 | (length(unique(hulldat$LatBin))==1)){
    #' 2. add a point to the south in the MIDDLE of the cluster.
    #' provided the shape of the survey area, seems like the most efficient option
    new_row = hulldat[hulldat$LatBin == min(hulldat$LatBin),]
    new_row$LatBin = min(hulldat$LatBin) - 0.2
    #' set middle of the area as Longitude of that row
    new_row$LonBin = mean(hulldat$LonBin)
    plot(hulldat$LonBin, hulldat$LatBin, pch=20)
    points(new_row$LonBin, new_row$LatBin, col="red", pch=20)
    hulldat = rbind(hulldat, new_row)
  }
  
  #' if ahull does function properly:
  if(i == 5){alpha=1.5}else{alpha=2} #' fix bug
  hulldat = hulldat %>% distinct() # remove duplicates
  alphahull_1 <- ashape(hulldat, alpha = alpha)
  plot(hulldat[,1:2], pch = 19, col = "darkseagreen")
  plot(alphahull_1, col = "magenta", add = TRUE)
  pts = hulldat[,1:2]%>%distinct()
  aoi <- ashape(pts, alpha = alpha)
  a <- data.frame(aoi$edges)[,c( 'x1', 'y1', 'x2', 'y2')]
  l <- st_linestring(matrix(as.numeric(a[1,]), ncol=2, byrow = T))
  for(k in 2:nrow(a)){
    l <- c(l, st_linestring(matrix(as.numeric(a[k,]), ncol=2, byrow = T)))
  }
  alphapoly <- st_sf(geom = st_sfc(l), crs = crs(x)) %>% 
    st_polygonize() %>% st_collection_extract()
  plot(hulldat[,1:2], pch = 19, col = "darkseagreen")
  plot(alphapoly,  add=T)
  xy3[xy3$Clust == i,]$geometry = alphapoly$geom
}

class(alphapoly$geom)

my_hulls = xy3 %>%
  arrange(-desc(LonBin))%>%
  group_by(geometry, nyears)%>%
  dplyr::summarise(Lon_centre=mean(LonBin),
                   Lat_centre=mean(LatBin))%>%
  ungroup()%>%
  dplyr::select(geometry, nyears, Lon_centre, Lat_centre)%>%
  distinct()%>%
  mutate(ID=1:length(unique(xy2$Clust)))


#' Plot final result
p1 = p+geom_sf(data=my_hulls, aes(geometry=geometry,  
                             fill=as.factor(ID), 
                             col=as.factor(ID)), alpha=0.3, lwd=1)+
  coord_sf(xlim = c(114, 133), ylim = c(-30, -36))+
  labs()+
  ggtitle(paste("dist =" , paste(max_dist2, "km", sep=" ")))+
  labs(fill="Cluster ID", col="Cluster ID")
p1

