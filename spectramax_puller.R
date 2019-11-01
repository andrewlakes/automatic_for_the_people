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
library(matrixStats)


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


####APPEND DATA####
endPoint = 187
plateArray = plateArray[,,1:endPoint]
timesOut = timesOut[1:endPoint,1]


#subtract background
#duplicate initial absorbance
plateArraybk = array(as.matrix(plateArray[,,1]), dim=c(8,12,endPoint))
plateArray = plateArray-plateArraybk

#divide by maximum absorbance found to normalize
plateArray = plateArray/max(plateArray)




#make into groups
#n=4, 3 groups


g1 = plateArray[,1:4,]
g2 = plateArray[,5:8,]
g3 = plateArray[,9:12,]


######Convert NA into the average in g1 since well messed up!#####
for (i in 1:length(plateArray[1,1,])){
    g1[3,4,i] = sum(g1[3,1,i], g1[3,2,i], g1[3,3,i])/3
}



gAve = array(NA, dim = c(8,3,length(timesOut)))
gErr = array(NA, dim = c(8,3,length(timesOut)))

for (i in 1:length(timesOut)){
  gAve[,1,i] = rowMeans(g1[,,i])
  gAve[,2,i] = rowMeans(g2[,,i])
  gAve[,3,i] = rowMeans(g3[,,i])
  
  gErr[,1,i] = rowSds(g1[,,i])
  gErr[,2,i] = rowSds(g2[,,i])
  gErr[,3,i] = rowSds(g3[,,i])
  
}


g1x = matrix(NA, nrow=8, ncol=1)
g2x = matrix(NA, nrow=8, ncol=1)
g3x = matrix(NA, nrow=8, ncol=1)

for (i in 1:7){
  g1x[i] = round(4*750/(4*(i^2)), digits = 0)
  g2x[i] = round(4*750/(4*(i^2)), digits = 0)
  g3x[i] = round(4*0.56/(4*(i^2)), digits = 3)
}

g1x[8] = 0
g2x[8] = 0
g3x[8] = 0

g1x = sapply(g1x, paste, " nCi/100uL")
g2x = sapply(g2x, paste, " nCi/100uL")
g3x = sapply(g3x, paste, " mg/mL")


pal = c("#1289d7", "#00008b", "#B20000", "#ffa500","#4C0CFC", "#0ACC1A", "#FF00A1", "#7F2A06")



plot1 = plot_ly() %>%
  add_trace(x = timesOut, y = gAve[1,1,], name = g1x[1], type = "scatter", mode = 'lines+markers', marker = list(color=pal[1]), line = list(color=pal[1])) %>%#, error_y = ~list(array = gErr[1,1,], color = pal[1])) %>%
  add_trace(x = timesOut, y = gAve[2,1,], name = g1x[2], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[2]), line = list(color=pal[2])) %>%#, error_y = ~list(array = gErr[2,1,], color = pal[2])) %>%
  add_trace(x = timesOut, y = gAve[3,1,], name = g1x[3], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[3]), line = list(color=pal[3])) %>%#, error_y = ~list(array = gErr[3,1,], color = pal[3])) %>%
  add_trace(x = timesOut, y = gAve[4,1,], name = g1x[4], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[4]), line = list(color=pal[4])) %>%#, error_y = ~list(array = gErr[4,1,], color = pal[4])) %>%
  add_trace(x = timesOut, y = gAve[5,1,], name = g1x[5], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[5]), line = list(color=pal[5])) %>%#, error_y = ~list(array = gErr[5,1,], color = pal[5])) %>%
  add_trace(x = timesOut, y = gAve[6,1,], name = g1x[6], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[6]), line = list(color=pal[6])) %>%#, error_y = ~list(array = gErr[6,1,], color = pal[6])) %>%
  add_trace(x = timesOut, y = gAve[7,1,], name = g1x[7], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[7]), line = list(color=pal[7])) %>%#, error_y = ~list(array = gErr[7,1,], color = pal[7])) %>%
  add_trace(x = timesOut, y = gAve[8,1,], name = g1x[8], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[8]), line = list(color=pal[8])) %>%#, error_y = ~list(array = gErr[8,1,], color = pal[8])) %>%
  layout(title = "DOTA-Ac-225",
         xaxis = list(title = "Time (hr)", showgrid = FALSE, range = c(0, 46)),
         yaxis = list(title = "Normalized Growth", showgrid = FALSE, range = c(0, 1)),
         font = list(size = 12))
  
plot1



plot2 = plot_ly() %>%
  add_trace(x = timesOut, y = gAve[1,2,], name = g2x[1], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[1]), line = list(color=pal[1])) %>%
  add_trace(x = timesOut, y = gAve[2,2,], name = g2x[2], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[2]), line = list(color=pal[2])) %>%
  add_trace(x = timesOut, y = gAve[3,2,], name = g2x[3], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[3]), line = list(color=pal[3])) %>%
  add_trace(x = timesOut, y = gAve[4,2,], name = g2x[4], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[4]), line = list(color=pal[4])) %>%
  add_trace(x = timesOut, y = gAve[5,2,], name = g2x[5], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[5]), line = list(color=pal[5])) %>%
  add_trace(x = timesOut, y = gAve[6,2,], name = g2x[6], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[6]), line = list(color=pal[6])) %>%
  add_trace(x = timesOut, y = gAve[7,2,], name = g2x[7], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[7]), line = list(color=pal[7])) %>%
  add_trace(x = timesOut, y = gAve[8,2,], name = g2x[8], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[8]), line = list(color=pal[8])) %>%
  layout(title = "IgG(#2)-DOTA-Ac-225",
         xaxis = list(title = "Time (hr)", showgrid = FALSE, range = c(0, 46)),
         yaxis = list(title = "Normalized Growth", showgrid = FALSE, range = c(0, 1)),
         font = list(size = 12),
         showlegend = TRUE)

plot2  
  
plot3 = plot_ly() %>%
  add_trace(x = timesOut, y = gAve[1,3,], name = g3x[1], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[1]), line = list(color=pal[1])) %>%
  add_trace(x = timesOut, y = gAve[2,3,], name = g3x[2], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[2]), line = list(color=pal[2])) %>%
  add_trace(x = timesOut, y = gAve[3,3,], name = g3x[3], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[3]), line = list(color=pal[3])) %>%
  add_trace(x = timesOut, y = gAve[4,3,], name = g3x[4], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[4]), line = list(color=pal[4])) %>%
  add_trace(x = timesOut, y = gAve[5,3,], name = g3x[5], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[5]), line = list(color=pal[5])) %>%
  add_trace(x = timesOut, y = gAve[6,3,], name = g3x[6], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[6]), line = list(color=pal[6])) %>%
  add_trace(x = timesOut, y = gAve[7,3,], name = g3x[7], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[7]), line = list(color=pal[7])) %>%
  add_trace(x = timesOut, y = gAve[8,3,], name = g3x[8], type = 'scatter', mode = 'lines+markers', marker = list(color=pal[8]), line = list(color=pal[8])) %>%
  layout(title = "IgG(#2)",
         xaxis = list(title = "Time (hr)", showgrid = FALSE, range = c(0, 46)),
         yaxis = list(title = "Normalized Growth", showgrid = FALSE, range = c(0, 1)),
         font = list(size = 12))

plot3

#subplot(plot1, plot2, plot3)



