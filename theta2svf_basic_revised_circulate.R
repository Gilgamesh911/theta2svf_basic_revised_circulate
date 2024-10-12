#conversion from Theta image to fish-eye image (equidistant and equal-area projection)
#binarization and calculatin of sky view factor
#theta2svf_basic, copyright 2018(c) Tsuyoshi Honjo
# 
library(jpeg)
library(png)

#following lines are changed if directory, file name, threshold value, rgb ratio are changed
setwd("C:/Users/1064112363/Desktop/eg") #directory of the Theta's image file
filelist <- list.files(getwd());#获取当前工作路径的文件列表
dir=paste("C:/Users/1064112363/Desktop/eg",filelist,sep="")#获取当前工作路径的文件列表的路径
N<-length(dir)#获取当前工作路径的文件数目
for(ii in 1:N)
{
  
  xth<-0.6 #threshold value in making binary images
  #rgb ratio
  r_ratio<-1
  g_ratio<-1
  b_ratio<-1
  
  #reading the data and initial process
  x<-readJPEG(filelist[ii],native=FALSE);#读取jpg 图形
  xdim<-nrow(x)
  x2<-array(0, dim=c(xdim+1, xdim+1, 3))
  x4<-x2
  xdim2<-xdim/2
  
  #transormation from cylindrical projection to equidistant and equal-area projection
  for(i in -xdim2:xdim2){
    for(j in -xdim2:xdim2){
      if(i>=0){ix<-xdim}
      else if (j>=0){ix<-xdim*2}
      else {ix<-0}
      ix<-round(atan(j/i)*xdim/pi)+ix
      iy<-sqrt(i^2+j^2)
      iy2<-round(xdim2*asin(iy/xdim2/sqrt(2))/pi*4) #equal-area
      iy<-round(iy) #eqidistance
      if (ix>=1 && ix<=xdim*2 && iy>=1 && iy<=xdim/2){
        x2[i+xdim2+1,j+xdim2+1,]<-x[iy,ix,]
      }
      if (ix>=1 && ix<=xdim*2 && iy2>=1 && iy2<=xdim/2){
        x4[i+xdim2+1,j+xdim2+1,]<-x[iy2,ix,]
      }
    }
  }
  
  #removing one line
  x2b<-x2[c(1:xdim2,(xdim2+2):(xdim2*2+1)),,]
  x2<-x2b[,c(1:xdim2,(xdim2+2):(xdim2*2+1)),]
  x4b<-x4[c(1:xdim2,(xdim2+2):(xdim2*2+1)),,]
  x4<-x4b[,c(1:xdim2,(xdim2+2):(xdim2*2+1)),]
  #saving fisheye images of equidistant (xxx_dist.png) and equal-area (xxx_area.png) 
  writePNG(x2,paste(ii,"_dist.png",sep=""))
  writePNG(x4,paste(ii,"_area.png",sep=""))
  
  #making equidistant binary image
  x<-(x2[,,1]*r_ratio+x2[,,2]*g_ratio+x2[,,3]*b_ratio)/(r_ratio+g_ratio+b_ratio)
  x[x>=xth]<-1
  x[x<xth]<-0
  #saving binary fisheye images of equidistant (xxx_dist_bi.png) 
  writePNG(x,paste(ii,"_dist_bi.png", sep=""))
  
  #making equal-area binary image
  x<-(x4[,,1]*r_ratio+x4[,,2]*g_ratio+x4[,,3]*b_ratio)/(r_ratio+g_ratio+b_ratio)
  x[x>=xth]<-1
  x[x<xth]<-0
  #saving binary fisheye images of equal-area (xxx_area_bi.png) 
  writePNG(x,paste(ii,"_area_bi.png", sep=""))
  
  #calculation of sky view factor from equal-area image
  (svf<-sum(x)/(pi*xdim2^2))
  write.table(svf,paste(ii,"_area_svf.txt", sep=""),col.names=F)
  
}
