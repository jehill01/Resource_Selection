library(adehabitatHR); library(stringr)
library(lme4); library(MuMIn); library(sp); library(rasterVis)
library(lubridate); library(rgdal); library(rgeos); library(broom)
library(raster); library(dplyr); library(rlist); library(pROC); library(ggplot2);
detach("package:ggsn", unload=TRUE)
data<-read.csv("dataNov.csv")
data$Seconds<-paste("0")
data$Timestamp<-paste(as.Date(data$Day, origin="2020/12/31"),"-",data$Hour,'-',data$Minutes,'-',data$Seconds)
data$Timestamp<-ymd_hms(data$Timestamp) #convert timestamp format
data<-with(data, data[(Timestamp <= "2021-07-31"),]) #only keeping the ones before this date
coordinates(data) <- c("coords.x1", "coords.x2") #defining coordinates from the input file
proj4string(data) <- CRS( "+proj=longlat +datum=WGS84" ) #setting the projection
data2<-spTransform(data, CRS('+proj=utm +zone=17N +datum=NAD83')) #reprojecting
animallist<-split(data2, data2$ID) #splitting by animal#

#plotting the points on a map of the site just to make sure things are looking right#
SITE<-readOGR("SITE.shp")
SITE2<-spTransform(SITE, CRS('+proj=utm +zone=17N +datum=NAD83'))
plot(SITE2)
colors<-as.list(c("red", "blue", "springgreen", "darkorange1", "pink", "steelblue1", "darkorchid1", 
                  "chartreuse3", "lightseagreen", "gold2", "blueviolet", "violetred4", "coral2",
                  "turquoise3", "tomato4", "darkgreen"))
for (i in 1:length(animallist)){
  points(animallist[[i]], col=colors[[i]], pch=20)
}

#running the kernels for every animal#

kernref<-list()
kernarea<-list()
kernelpoly95<-list()
kernelpoly25<-list()
polydiff<-list()
for (i in seq_along(animallist)) {
  kernref [[i]]<- kernelUD(animallist[[i]], h = "href")
  kernarea[[i]]<-kernel.area(kernref[[i]], percent = 95,
              unin = c("m"),
              unout = c("km2"), standardize = FALSE) #gives the area, don't really need for this
  kernelpoly95[[i]]<-getverticeshr(kernref[[i]], percent= 95, unin = c("m"), unout = c("km2")) #95% kernel
  kernelpoly25[[i]]<-getverticeshr(kernref[[i]], percent= 25, unin = c("m"), unout = c("km2")) #25% kernel
  polydiff[[i]]<-gDifference(kernelpoly95[[i]], kernelpoly25[[i]]) #subtracting 25% from the 95% 
  polydiff[[i]]$ID<-paste(animallist[[i]]$ID[1])
  polydiff[[i]]$Sex<-paste(animallist[[i]]$Sex[1])
  polydiff[[i]]$Habitat<-paste(animallist[[i]]$Habitat[1])
}



available<-list()
used<-list()
RSFpts<-list()

for (i in seq_along(polydiff)) {
  available[[i]]<-spsample(polydiff[[i]], n=75, type="random")@coords #sample from the 95-25% kernel for available
  used[[i]]<-spsample(kernelpoly25[[i]], n=75, type="random")@coords  #sample from the 25% for used
  RSFpts[[i]] <- data.frame(coords.x1 = c(used[[i]][, 1], available[[i]][, 1]), 
                            coords.x2 = c(used[[i]][, 2], available[[i]][, 2]), 
                            Used = c(rep(TRUE, nrow(used[[i]])), rep(FALSE, nrow(available[[i]]))))
  RSFpts[[i]]$ID<-paste(polydiff[[i]]$ID) #appending all the animal info
  RSFpts[[i]]$Sex<-paste(polydiff[[i]]$Sex)
  RSFpts[[i]]$Habitat<-paste(polydiff[[i]]$Habitat)
}

lcnew<-raster("Reclass_nlcd11.tif") #reading in our land cover raster
merger<-do.call(rbind, RSFpts) 
combined<-cbind(merger$coords.x1, merger$coords.x2)

out<-list()
land<-list()
merge<-list()
mergedf<-list()

for (i in 1:length(RSFpts)) {
  land[[i]]<-as.data.frame(extract(lcnew, RSFpts[[i]][,1:2])) #extracting the land cover, not really using for this
  merge[[i]]<-list.merge(as.data.frame(RSFpts[[i]]), as.data.frame(land[[i]]))
  mergedf[[i]]<-as.data.frame(merge[[i]])
}

#This gets the summary stats of used vs available habitats
all<-do.call(rbind, mergedf)
names(all)[names(all) == "extract.lcnew..RSFpts..i.....1.2.."] <- "LandCover"
all<- all %>% mutate(Cover =
                               case_when(LandCover==0~ "OpenWater",
                                 LandCover==1~ "Developed",
                                 LandCover==2~ "Barren",
                                 LandCover==3~ "Deciduous",
                                 LandCover==4~ "Evergreen",
                                 LandCover==5~ "MixedForest",
                                 LandCover==6~ "Shrub",
                                 LandCover==7~"Grass",
                                 LandCover==8~ "Agriculture",
                                 LandCover==9~ "WoodyWetlands",
                                 LandCover==10~ "HerbaceousWetlands"))

LCSTATS<-all %>% 
  group_by(Used) %>% count(Cover)  
LCSTATS <- arrange(LCSTATS, Cover, Used)
LCSTATS%>% print(n=50)
#Distance from every point to every land cover type
for (i in 1:nrow(combined)) {#calculating distance from every point to every lc type
  #this takes a while to run
  d <- distanceFromPoints(lcnew, combined[i,,drop=F])
  out[[i]] <- zonal(d, lcnew, min)[,2]
}

a<-do.call(rbind, out)
b<-data.frame(a)
c<-bind_rows(RSFpts)
RSF<-bind_cols(b,c)
colnames(RSF)<-c("OpenWater", "Developed", "Barren", "Deciduous", "Evergreen", "Mixedforest", "Shrub", "Grass",
               "Agriculture", "Woodywetlands", "Herbwetlands", "coords.x", "coords.y", "Used", "ID", "Sex", "Habitat")

RSF<-readRDS("distance.RDS") #RSF saved here, reading it in when resuming

RSF<-subset(RSF, select=-c(OpenWater,Agriculture, Herbwetlands, Mixedforest, Shrub, Grass, Barren)) #Dropping these because they are sparse

#doing a correlation matrix to see which land cover types are correlated
corrmatrix<-RSF[ , -which(names(RSF) %in% c("coords.x","coords.y", "Used", "ID", "Sex", "Habitat"))] #Remove everything but lc distances
cor(corrmatrix, method = "pearson", use = "complete.obs") #run the correlation
#nothing >0.70 cutoff so using all four lc types
fullmodel<-glm(Used ~scale(Woodywetlands)+ scale(Deciduous)+scale(Evergreen)+scale(Developed), data=RSF, family="binomial", na.action=na.pass)
dredge(fullmodel)
#based on dredge, full model was most supported so keeping all lc types

#This table provides mean distance to each land cover type
STATS<-RSF %>% 
  group_by(Used) %>% 
  summarise(Meandev=mean(Developed),
            Meandecid=mean(Deciduous),
            Meanever=mean(Evergreen),
            Meanwoody=mean(Woodywetlands)
   )



auc(Used~predict(fullmodel), data=RSF, plot=TRUE) #calculating area under curve


#Upload the lc rasters, created in other R file#
woody<-raster("Woodywetlands.tif")
ever<-raster("Evergreen.tif")
Deciduous<-raster("Deciduous.tif")
Developed<-raster("Developed.tif")

#stack all the landcovers 
landraster<-stack(c(woody, ever, Deciduous, Developed))

#generate prediction raster using raster stack and model
predictmap<-predict(landraster, fullmodel, fun=predict, na.rm=TRUE, type="response")
plot(predictmap)
##Figures
levelplot(predictmap, margin=FALSE, main='Predicted Probability of Use')

RSF2<-RSF %>% mutate(Treatment =case_when(Used=="TRUE" ~ "Used",
                                          Used=="FALSE" ~ "Not used"))

ggplot(RSF2,aes(x=Evergreen))+geom_histogram(binwidth = 20)+facet_grid(~Treatment)+theme_bw()+
  xlab("Distance to evergreen land cover")+ylab("Cell count")
ggplot(RSF2,aes(x=Woodywetlands))+geom_histogram(binwidth=20)+facet_grid(~Treatment)+theme_bw()+
  xlab("Distance to woody wetland land cover")+ylab("Cell count")

#Example plot of used/available
plot(polydiff[[1]])
scalebar(1000, type="bar", divs=4, below="Meters", adj=c(0.5, -1.0))
points(available[[1]], add=TRUE, col="blue", pch=19)
points(used[[1]], add=TRUE, col="red", pch=19)

