library(reshape2)
library(plyr)
library(dplyr)
library(drc)
library(drfit)
library(ggplot2)
library(deSolve)
library(emdbook)
library(stats)
library(plotly)
library(cowplot)
library(gridExtra)
library(abind)
library(RColorBrewer)
library(tidyr)
library(GenKern)
library(xlsx)
library(ggpubr)
library(grid)
library(lubridate)


dftabs = read.delim('IVK_plague_46hr_final.pda.txt', header = FALSE)

dftabs = dftabs[-c(1,2),]

times = as.data.frame(dftabs[2:length(dftabs[,1]),1])
times = as.data.frame(times[-c(nrow(times)-1,nrow(times)),])
times = as.data.frame(times[seq(1,nrow(times),9),1])



timesOut = as.data.frame(matrix(NA, nrow = length(times[,1]), ncol = 1))


for (i in 1:length(times[,1])){
  if (nchar(toString(times[i,1])) < 7) {timesOut[i,1] = period_to_seconds(ms(toString(times[i,1])))}
  else {timesOut[i,1] = period_to_seconds(hms(toString(times[i,1])))}
}

#convert from seconds to hours
timesOut = timesOut/60/60
colnames(timesOut) = "Time (hr)"



plate = as.matrix(dftabs[2:length(dftabs[,1]),3:14])
colnames(plate) = 1:12


#convert plate into an array for each time

#96 well plate
plateArray = array(NA, dim = c(8,12,length(timesOut[,1])))

for (i in 1:length(timesOut[,1])){
  plateArray[,,i] = as.numeric(plate[seq(1+9*(i-1),8+9*(i-1),1),])
}


plateArray = plateArray[,,1:187]
timesOut = timesOut[1:187,1]


#subtract background
plateArray = plateArray - 0.056

#divide by maximum absorbance found to normalize
plateArray = plateArray/max(plateArray)

#make into groups
#n=4, 3 groups


g1 = plateArray[,1:4,]
g2 = plateArray[,5:8,]
g3 = plateArray[,9:12,]

gAve = array(NA, dim = c(3,8,length(timesOut)))
gErr = array(NA, dim = c(3,8,length(timesOut)))

for (i in 1:length(timesOut)){
  gAve[,1,i] = rowMeans(g1[,,i])
  gAve[,2,i] = rowMeans(g2[,,i])
  gAve[,3,i] = rowMeans(g3[,,i])
  
  gErr[,1,i] = rowSds(g1[,,i])
  gErr[,2,i] = rowSds(g2[,,i])
  gErr[,3,i] = rowSds(g3[,,i])
  
}










mAverage225tras = melt(Average225tras, id="Days")
colnames(mAverage225tras) = c("times", "Organs", "values")

plot225tras = ggplot()+ 
  geom_line(data=mAverage225tras, aes(x=times, y=values, color=Organs), size=1, alpha=1)+
  geom_ribbon(data=mAverage225errortras, aes(x=times, ymin=valuesminus,  ymax=valuesplus, fill = Organs), alpha = 0.1)+
  ggtitle("DOTA-Trastuzumab")+
  
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))+ 
  
  scale_x_log10()+#breaks=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000))+
  annotation_logticks(base = 10, sides = "b", scaled = TRUE,
                      short = unit(0.1, "cm"), mid = unit(0.2, "cm"), long = unit(0.3, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL)+
  
  scale_y_continuous()+#limits = c(min(plot225scale),max(plot225scale)), breaks=plot225scale)+#breaks=c(lseq(0.000001,100,9)))+
  theme_bw() +
  theme(legend.position="none", plot.margin = unit(margins, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y=element_blank(), x=element_blank(), color="Organs")+
  theme(text = element_text(size=18, face = "bold"),
        axis.text.y=element_text(colour="black"),
        axis.text.x=element_text(colour="black"),
        plot.title = element_text(hjust = 0.5, size=18))
#+
#guides(shape=guide_legend(override.aes = list(size=3)))
