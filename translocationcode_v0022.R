rm(list=ls(all=TRUE))
library(ggeffects); library(sjPlot)
library(segclust2d); library(move); library(openair)
library(dplyr); library(lme4); library(MuMIn)
library(ggplot2);library(emmeans)
set.seed(34)


##Read in control animals###
controls<-read.csv("precontrolcollar.csv")
unique(controls$collarID)
controls$Timestamp<-as.POSIXct(controls$Timestamp, tz="EST", 
	format="%m/%d/%Y %H:%M:%S") #convert time to POSIX#
###Sampling set period of time (to match transient/resident time periods)###
controls$date<-controls$Timestamp 
	#selectByDate function requires column named 'date' (lowercase!!) #
controlwint<-selectByDate(controls, 
	month = c("December","January","February", 
			"November","October","March"))
newdata <- controlwint[order(controlwint$collarID, controlwint$Timestamp),]
	 #ordering the dataframe#
contlist<-split(newdata, newdata$collarID) 
	#split into list by collarID#
monthlist<-list() #make a list for the loop#
names(contlist)
for (i in 1:length(contlist)) {
  monthlist[[i]]<-cut(contlist[[i]]$Timestamp, breaks="45 days")
  
} #break up each animal into 45 days chunks"


timebreaks<-data.frame(unlist(monthlist)) #unlist and put into dataframe#
timebreaks$rown <- 1:nrow(timebreaks) #add column with row number#
newdata$rown<-1:nrow(newdata) #add column with row number#
dataset<-merge.data.frame(newdata, timebreaks, 
	by.y="rown", by.x="rown") 
	#merge dataset with time breaks based 
	#on row numbers (won't work if both not in same order)#
dataset$ID_Timebreak<-paste(dataset$collarID, 
	dataset$`unlist.monthlist.`) 
	#add a column with every time break for every animal#
duration<-dataset %>% group_by(ID_Timebreak) %>% 
	summarise(firstvis=min(Timestamp), 
	lastvis=max(Timestamp), duration=lastvis-firstvis) 
	#tibble with the duration of every time break (not all div by 45)#
complete<-merge.data.frame(duration, dataset, 
	by.y="ID_Timebreak", by.x="ID_Timebreak") 
	#merge duration with the dataframe#
complete %>% filter(duration > 44) #filter out ones less than 45#
## complete[5000:6000,]
by_cyl <- complete %>% group_by(collarID) 
	#grouping by ID for time selection#
r77<-sample_n(by_cyl, 1, replace = FALSE) 
	#sampling one time break for every animal#
controlfinal<-complete %>%
  filter(ID_Timebreak %in% r77$ID_Timebreak) 
	#filter the dataframe based on the breaks selected#
## controlfinal[5000:5500,] #dataframe merging can get wonky 
	# so spot check to make sure everything's lined up"
coordinates(controlfinal)<-c("Easting","Northing")
proj4string(controlfinal)<-CRS("+proj=longlat +datum=WGS84")
control2<-as.data.frame(spTransform(controlfinal, 
	CRS('+proj=utm +zone=17N +datum=NAD83'))) #convert to UTM#
control2$Timestamp<-as.POSIXct(control2$Timestamp, 
	tz="EST", format="%m/%d/%Y %H:%M:%S")

###Read in translocated and perform segmentation###

moved<-read.csv("collarstranslocation.csv")
moved<- moved %>% distinct() #removing duplicates#
coordinates(moved)<-c("Easting","Northing")
proj4string(moved)<-CRS("+proj=longlat +datum=WGS84")
moveddf<-as.data.frame(spTransform(moved, 
	CRS('+proj=utm +zone=17N +datum=NAD83')))
moveddf$Timestamp<-as.POSIXct(moveddf$Timestamp, 
	tz="EST", format="%m/%d/%Y %H:%M:%S")
movedlist<-split(moveddf, moveddf$ID)
unique(moved$collarID)
finalmoved<-movedlist

finalmoved<-movedlist[names(movedlist) %in% 
                        c("40996 PostMove", "41011 PostMove", "41002 PostMove") == FALSE] 
		#dropping data deficient ones#

shift.list<-list() #create list for the loop#
shift.info<-list() #create list for the loop#
for (i in seq_along(finalmoved)) {
  shift.list[[i]]<-segmentation(finalmoved[[i]],lmin=25,  
		Kmax = 7, seg.var = c("Easting","Northing"), 
		scale.variable=FALSE) 
  shift.info[[i]]<-augment(shift.list[[i]])
  
}


segresults <- do.call("rbind", shift.info) 
	#merge list elements back into a single dataframe#
segresults$ID_state<-paste(segresults$collarID, 
	segresults$state) #create column with ID for 
				#every state for every animal#

#back to the original task#
trans_res<-segresults %>% filter(state <= 2 ) #only keeping 1st two states#
transkeep=trans_res[c("collarID","SEX","state","Treatment", 
	"PrePost","Timestamp","Easting","Northing","ID_state")]  
	#keep only the needed columns#

control2$state<-0 #Assigning controls state of 0 to 
			#match state of 1 for trans and 2 for res#
control2$ID_state<-paste(control2$collarID, control2$state)

controlkeep=control2[c("collarID","SEX","state","Treatment", 
	"PrePost","Timestamp","Easting","Northing","ID_state")] 
	#make same columns as trans for binding"

finaldata<-rbind(controlkeep, transkeep) #bind the two dataframes"

finaldata$Group[finaldata$state==0]<-"Control" 
	#Making new column and assigning control/trans/res 
	#based on state number#
finaldata$Group[finaldata$state==1]<-"Transient" 
finaldata$Group[finaldata$state==2]<-"Resident"
## head(finaldata)
final2<-finaldata[with(finaldata, order(collarID, Timestamp)),] 
	#need ascending timestamps for every individual"
final2<-final2[!(final2$collarID=="41015"),]
final2<-final2[!(final2$collarID=="41002"),]
final2<-final2[!(final2$collarID=="40996"),]
unique(final2$collarID)
###Home range construction###
raccoons<-move(x=final2$Easting, y=final2$Northing, 
	animal=final2$ID_state, data=final2, time=final2$Timestamp)
datalist<-split(raccoons, final2$ID_state) ## names(datalist)

modlist95<-list()
udlist95<-list()
plotlist95<-list()
arealist95<-list()
modlist60<-list()
udlist60<-list()
plotlist60<-list()
arealist60<-list()

for (i in 1:length(datalist) ){ ## i=names(datalist)[2]
  modlist95[[i]]<-brownian.bridge.dyn(datalist[[i]], location.error=7, raster=50, ext=4, margin=5, window=11)
  flush.console()
  udlist95[[i]]<-getVolumeUD(modlist95[[i]])
  plotlist95[[i]]<-udlist95[[i]]<=0.95
  arealist95[[i]]<-sum(values(plotlist95[[i]]))
  
 
}


getVolumeUD(modlist95[[1]])

brownian.bridge.dyn(datalist[[12]],location.error=7, raster=50, ext=c(5,5,6,6), margin=5, window=13)

plot(plotlist95[[41]], ext=c( 440000, 460000, 3660000, 3690000))
points(datalist[[41]])
data.frame(datalist[[48]]$Timestamp)

#Ouput of arealist is number of raster cells, so
	# multiply area by the size of the raster cells then 
	# divide by 1e06 to get sq km#
Area<-lapply(arealist95, function(x) (x*2500)/1000000) 
Area2<-as.data.frame(unlist(Area)) #get the areas into a df#
Individual<-as.data.frame(names(datalist))
AreaDF<-merge.data.frame(data.frame(cbind(Area2, Individual)), 
	finaldata, by.x="names.datalist.", by.y="ID_state" )

modlist60<-list()
udlist60<-list()
plotlist60<-list()
arealist60<-list()

for (i in names(datalist) ){ ## i=1
  modlist60[[i]]<-brownian.bridge.dyn(datalist[[i]], 
	location.error=7, raster=50, ext=4, margin=5, window=13,verbose=F)
  flush.console()
  udlist60[[i]]<-getVolumeUD(modlist60[[i]])
  plotlist60[[i]]<-udlist60[[i]]<=0.60
  arealist60[[i]]<-sum(values(plotlist60[[i]]))
}

Area60<-lapply(arealist60, function(x) (x*2500)/1000000) 
Area260<-as.data.frame(unlist(Area60))
AreaDF2<-merge.data.frame(data.frame(cbind(Area260, Individual)), 
	AreaDF, by.x="names.datalist.", by.y="names.datalist." )
AreaDF2


AreaDF3<-AreaDF2
names(AreaDF3)[names(AreaDF3) == "unlist.Area."] <- "UD95"
names(AreaDF3)[names(AreaDF3) == "unlist.Area60."] <- "UD60"
names(AreaDF3)[names(AreaDF3) == "names.datalist."] <- "Individual"
AreaDF3<-AreaDF3[!(AreaDF3$Individual=="40997 2"),]
AreaDF3<-AreaDF3[!(AreaDF3$Individual=="41013 2"),]
AreaDF3=AreaDF3[order(AreaDF3$Individual,AreaDF3$Timestamp),]

AreaDF3<- AreaDF3 %>% mutate(Start =
	case_when(
		Treatment=="H2L"  ~ "High",
		Treatment=="H2H"  ~ "High", 
		Treatment=="L2L"  ~ "Low",
		Treatment=="L2H"  ~ "Low",
		collarID=="40993"~ "Low",
		collarID=="41001"~ "Low",
		collarID=="41010"~ "Low",
		collarID=="41012"~ "High",
		collarID=="41015"~ "High",
		collarID=="50000"~ "Low",
		collarID=="70000"~ "Low"))


AreaDF3<- AreaDF3 %>% mutate(End =
	case_when(
		Treatment=="H2L"  ~ "Low",
		Treatment=="H2H"  ~ "Low", 
		Treatment=="L2L"  ~ "Low",
		Treatment=="L2H"  ~ "High",
		collarID=="40993"~ "Low",
		collarID=="41001"~ "Low",
		collarID=="41010"~ "Low",
		collarID=="41012"~ "High",
		collarID=="41015"~ "High",
		collarID=="50000"~ "Low",
		collarID=="70000"~ "Low"))

saveRDS(AreaDF3, "AreaDF3a.rds")
AreaDF3<-readRDS("AreaDF3a.rds")
AreaInput2<-subset(AreaDF3, !duplicated(Individual)) 
	#keep one row per ID_state to tidy things up#

AreaDF3

unique(AreaDF3$Individual)

AreaInput2 %>%  group_by(state) %>% 
  tally() %>% print(n=27)
nrow(AreaInput2)

AreaDF3<-AreaDF3[(AreaDF3$PrePost=="PreMove"),]
AreaDF3<-AreaDF3[(AreaDF3$Start=="High"),]
duration<-AreaDF3 %>% group_by(Individual) %>% 
  summarise(firstvis=min(Timestamp), 
            lastvis=max(Timestamp), duration=lastvis-firstvis) 
sum(duration$duration)

##-------------------------------------------------------------------------------------
while(.Device != "null device") dev.off()
AreaInput<-AreaDF3
COLS=c( "collarID" ,"SEX","Start" ,"End")
idVec= unique(AreaInput$Individual); 
idVec=idVec[substr(idVec,nchar(idVec),nchar(idVec)) != "0"]
idVec=as.character(sapply(idVec,function(x) strsplit(x," ")[[1]][1]))
idVec=idVec[duplicated(idVec)]

distCentroid.DF=NULL
for(id in idVec){## id="80000"
  xy1= aggregate(extent(modlist60 [[ paste(id,"1") ]])[],
	by=list(c(1,1,2,2)),mean)$x
  xy2= aggregate(extent(modlist60 [[ paste(id,"2") ]])[],
	by=list(c(1,1,2,2)),mean)$x

  Dist = pointDistance(xy1,xy2,lonlat=F)/1000
  ROW=which(AreaInput$collarID==id)[1]
  DF=data.frame(AreaInput[ROW,COLS],distTR.km=Dist); rownames(DF)=NULL
  distCentroid.DF=rbind(distCentroid.DF,DF); 
  rm(xy1,xy2,Dist,DF)
}


##-------------------------------------------------------------------------------------
head(AreaDF3)
COLS=c("Individual","collarID","SEX","state","Treatment","Group",
	"Start" ,"End")
distHourly.DF=NULL
for(id in unique(AreaDF3$Individual)){ ## id=unique(AreaDF3$Individual)[1]
  ROW=which(AreaDF3$Individual==id)
  DF=AreaDF3[ROW,c("Timestamp","Easting","Northing")]
  DF$hours= c(999,diff(DF$Timestamp))
  R2= which( DF$hours <= 4 & DF$hours > 0); 
  Dist=sapply(1:length(R2),function(x) ## x=1
	pointDistance(as.numeric(DF[ R2[x]-1, c("Easting","Northing") ]),
		as.numeric(DF[ R2[x], c("Easting","Northing") ]),lonlat=F)/1000)
  DF=data.frame(AreaDF3[rep(ROW[1],length(Dist)),COLS],
	hours=DF$hours[R2],dist.km=Dist)
  DF$kmph=DF$dist.km/DF$hours
  row.names(DF)=NULL
  distHourly.DF=rbind(distHourly.DF,DF); 
  rm(Dist,DF,ROW,R2)
}

##-------------------------------------------------------------------------------------

COLS=c("Individual","collarID","SEX","state","Treatment","Group",
	"Start" ,"End")
dist12h.night.DF=dist12h.day.DF=dist24h.den.DF=NULL
for(id in unique(AreaDF3$Individual)){ ## id=unique(AreaDF3$Individual)[1]
  ROW=which(AreaDF3$Individual==id)
  DF=AreaDF3[ROW,c("Timestamp","Easting","Northing")]
  DF$HR=format(DF$Timestamp,"%H")
  DF=DF[ DF$HR %in% c("06","18"),]
  x=diff(DF$Timestamp); units(x)<-"hours"
  DF$hours= c(999,x)
  R2= which( DF$hours == 12); R2=R2[ DF$HR[R2] =="06"]
  ## DF[R2,] ; DF[R2-1,]
  Dist=NULL
  if(length(R2)>0){
    Dist=sapply(1:length(R2),function(x) ## x=1
	pointDistance(as.numeric(DF[ R2[x]-1, c("Easting","Northing") ]),
		as.numeric(DF[ R2[x], c("Easting","Northing") ]),lonlat=F)/1000)
    DF=data.frame(AreaDF3[rep(ROW[1],length(Dist)),COLS],
	hours=DF$hours[R2],dist.km=Dist)
    row.names(DF)=NULL
    dist12h.night.DF=rbind(dist12h.night.DF,DF); 
  }
  rm(Dist,DF,ROW,R2)

  ##---------------------------------------
  ROW=which(AreaDF3$Individual==id)
  DF=AreaDF3[ROW,c("Timestamp","Easting","Northing")]
  DF$HR=format(DF$Timestamp,"%H")
  DF=DF[ DF$HR %in% c("06","18"),]
  x=diff(DF$Timestamp); units(x)<-"hours"
  DF$hours= c(999,x)
  R2= which( DF$hours == 12); R2=R2[ DF$HR[R2] =="18"]
  ## DF[R2,] ; DF[R2-1,]
  Dist=NULL
  if(length(R2)>0){
    Dist=sapply(1:length(R2),function(x) ## x=1
	pointDistance(as.numeric(DF[ R2[x]-1, c("Easting","Northing") ]),
		as.numeric(DF[ R2[x], c("Easting","Northing") ]),lonlat=F)/1000)
    DF=data.frame(AreaDF3[rep(ROW[1],length(Dist)),COLS],
	hours=DF$hours[R2],dist.km=Dist)
    row.names(DF)=NULL
    dist12h.day.DF=rbind(dist12h.day.DF,DF); 
  }
  rm(Dist,DF,ROW,R2)

  ##-----------------------------------------
  ROW=which(AreaDF3$Individual==id)
  DF=AreaDF3[ROW,c("Timestamp","Easting","Northing")]
  DF$HR=format(DF$Timestamp,"%H")
  DF=DF[ DF$HR %in% c("12"),]
  x=diff(DF$Timestamp); units(x)<-"hours"
  DF$hours= c(999,x)
  R2= which( DF$hours == 24); 
  ## DF[R2,] ; DF[R2-1,]
  Dist=NULL
  if(length(R2)>0){
    Dist=sapply(1:length(R2),function(x) ## x=1
	pointDistance(as.numeric(DF[ R2[x]-1, c("Easting","Northing") ]),
		as.numeric(DF[ R2[x], c("Easting","Northing") ]),lonlat=F)/1000)
    DF=data.frame(AreaDF3[rep(ROW[1],length(Dist)),COLS],
	hours=DF$hours[R2],dist.km=Dist)
    row.names(DF)=NULL
    dist24h.den.DF=rbind(dist24h.den.DF,DF); 
  }
  rm(Dist,DF,ROW,R2)

}
dist24h.den.DF
saveRDS(AreaInput, paste0("AreaInput3.rds"))
saveRDS(AreaDF3, paste0("AreaDF3.rds"))
saveRDS(final2, paste0("final2.rds"))
saveRDS(distCentroid.DF,paste0("distCentroid.DF.rds"))
saveRDS(distHourly.DF,paste0("distHourly.DF.rds"))
saveRDS(dist12h.night.DF,paste0("dist12h.night.DF.rds"))
saveRDS(dist12h.day.DF,paste0("dist12h.day.DF.rds"))
saveRDS(dist24h.den.DF,paste0("dist24h.den.DF.rds"))

AreaInput
distHourly.DF<-readRDS("distHourly.DF.rds")
head(distHourly.DF)
head(dist24h.den.DF)
unique(dist12h.night.DF$collarID)
unique(AreaInput$collarID)

AreaInput2$Group <- factor(AreaInput2$Group, levels=c("Control", "Transient", "Resident"))


UD95full<-lm(log(UD95)~(Group+SEX+Start+End)^2, data=AreaInput2, na.action = na.pass)
UD60full<-lm(log(UD60)~(Group+SEX+Start+End)^2, data=AreaInput2, na.action = na.pass)
nightdist<-lmer(log(dist.km)~(Group+SEX+Start+End)^2+(1|collarID),data=dist12h.night.DF, na.action = na.pass, REML=F)
daydist<-lmer(log(dist.km)~(Group+SEX+Start+End)^2+(1|collarID),data=dist12h.day.DF, na.action = na.pass, REML=F)
hourdist<-lmer(log(kmph)~(Group+SEX+Start+End)^2+(1|collarID),data=distHourly.DF, na.action = na.pass, REML=F)
dendist<-lmer(log(dist.km)~(Group+SEX+Start+End)^2+(1|collarID),data=dist24h.den.DF, na.action = na.pass, REML=F)
dredge(UD95full)
dredge(UD60full)
dredge(hourdist)
dredge(nightdist)
dredge(daydist)

dredge(dendist)



UD95fit<-lm(log(UD95)~Group+SEX+Start, data=AreaInput2, na.action = na.pass)
ggpredict(UD95fit, ci.lvl=0.95, terms=c("Start"), condition = c(SEX="M", Group="Resident"))
ggpredict(UD95fit, ci.lvl=0.95, terms=c("SEX"), condition = c(Start="Low", Group="Resident"))
ggpredict(UD95fit, ci.lvl=0.95, terms=c("Group"))
emmeans(UD95fit, pairwise~Group, type="response")

UD60fit<-lm(log(UD60)~Group+SEX+Start+Group:Start, data=AreaInput2, na.action = na.pass)
ggpredict(UD60fit, ci.lvl=0.95, terms=c("SEX"), condition = c(Start="Low", Group="Resident"))
ggpredict(UD60fit, ci.lvl=0.95, terms=c("Group"), condition = c(SEX="M", Start="Low"))
ggpredict(UD60fit, ci.lvl=0.95, terms=c("Group"), condition = c(SEX="M", Start="High"))

emmeans(UD60fit, pairwise~Start|Group, type="response")
emmip(UD60fit, ~Start|Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")

nightfit<-lmer(log(dist.km)~Group+SEX+Start+End+Group:SEX+Group:Start+(1|collarID),data=dist12h.night.DF, na.action = na.pass, REML=F)
emmeans(nightfit, pairwise~Group|Start, type="response")
emmip(nightfit, ~Group|Start, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")
emmeans(nightfit, pairwise~SEX|Group, type="response")
emmip(nightfit, ~SEX|Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")





STAT<-AreaInput2 %>% 
  group_by(SEX, Group, Start) %>% 
  summarise(Mean95=mean(UD95),
            sd95=sd(UD95),
            mean60=mean(UD60),
            sd60=sd(UD60),
             )
STAT

STAT<-AreaInput2 %>% 
  summarise(Mean95=mean(UD95),
            sd95=sd(UD95),
            mean60=mean(UD60),
            sd60=sd(UD60),
  )
STAT
unique(AreaInput2$collarID)

UD60lm<-lm(log(UD60)~(Group+SEX+Start+End)^2, data=AreaInput2, na.action = na.pass)
dredge(UD60lm)
summary(UD95fit)
predict(UD95fit, data.frame(Group="Transient", SEX="M", Start="Low"), se.fit=TRUE, type=c('response'))
predict(UD95fit, data.frame(Group="Transient", SEX="M", Start="High"), se.fit=TRUE, type=c('response'))


emmip(fit95, ~Group, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")
emmeans(UD95fit, pairwise~Group, type="response")

emmip(fit95, ~Group|SEX, CIs=TRUE, CIarg=aes(lwd=1, color="black", linetype=1), type="response")+
  facet_wrap(~SEX)+
  theme(strip.background=element_rect(color="black", fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", 
        axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-3),
        axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.9,0.9), "cm"))+xlab("State")+ylab("95% UD (kilometers squared)")


plot_model(fit95, type="pred", terms = c("Group", "SEX", "Start"))+theme(strip.background=element_rect(color="black", fill="lightblue"), panel.border = element_rect(colour='black', fill='NA'), 
                                                                         panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
                                                                         legend.key = element_rect(fill='white'), legend.position="bottom", 
                                                                         axis.text=element_text(size=12), axis.title.x = element_text(size=14, vjust=-3),
                                                                         axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=12, face="bold"))+
  theme(plot.margin=unit(c(0.5,0.5,0.9,0.9), "cm"))+xlab("State")+ylab("95% UD (kilometers squared)")
f95<-lmer(log(UD95)~SEX+(1|collarID),data=AreaInput2, na.action = na.pass, REML=F)

ggpredict(UD95fit, ci.lvl=0.95, terms=c("Start"), condition = c(SEX="M", Group="Resident"))
ggpredict(fit95, ci.lvl=0.95, terms=c("Group", "SEX [M,F]"), condition = c(Start="Low"))
ggpredict(fit60, ci.lvl=0.95, terms=c("Group", "SEX", "Start"))

emmeans(fit60, pairwise~Group|Start, type="response")



 

