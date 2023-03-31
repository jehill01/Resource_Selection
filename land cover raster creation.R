
d4<-raster("Reclass_nlcd11.tif")
#reclassify each raster so the lc of interest is NA
evermatrix<-cbind(c(0:10),c(NA,NA,NA,NA,4,NA,NA,NA,NA,NA,NA))
everrecl<-reclassify(d4, rcl=evermatrix)
plot(everrecl)
everdist<-distance(everrecl)
names(everdist)<-"Evergreen" #name of layer must match name in data file

decidmatrix<-cbind(c(0:10),c(NA,NA,NA,3,NA,NA,NA,NA,NA,NA,NA))
decidrecl<-reclassify(d4, rcl=decidmatrix)
plot(decidrecl)
deciddist<-distance(decidrecl)
names(deciddist)<-"Deciduous" #name of layer must match name in data file

devmatrix<-cbind(c(0:10),c(NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA))
devrecl<-reclassify(d4, rcl=devmatrix)
plot(devrecl)
devdist<-distance(devrecl)
names(devdist)<-"Developed" #name of layer must match name in data file