# #betaSPM = read.csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/SPMbeta_num_obs_50sigma_1dist_noise_normali_1.csv")[,-1]
# setwd("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data")
# num_obs = 50
# sigma = 2
# i = 2
# dist_noise = "normal"
# voxels = read.csv("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/cube_data/voxels_num_obs_100sigma_6dist_noise_normali_4.csv")
# voxels = as.matrix(voxels)
# 
# true_beta = read.csv(paste0("true_beta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# true_beta = as.matrix(true_beta)
# 
# SPMbeta = read.csv(paste0("SPMbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# SPMbeta = as.matrix(t(SPMbeta[,-1]))
# 
# MUAbeta = read.csv(paste0("MUAbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# MUAbeta = as.matrix(t(MUAbeta[,-1]))
# 
# BSTbeta = read.csv(paste0("BSTbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# BSTbeta = as.matrix(t(BSTbeta[,-1]))
# 
# IOSMEbeta = read.csv(paste0("est_beta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# IOSMEbeta = as.matrix(IOSMEbeta)
# 
# IRRNNbeta = read.csv(paste0("IRRNNbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))[,-1]
# IRRNNbeta = t(as.matrix(IRRNNbeta))
# 
# STOREbeta = read.csv(paste0("STOREbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))
# STOREbeta = as.matrix(STOREbeta)
# 
# SVCMbeta = read.csv(paste0("SVCMbeta_", "num_obs_", num_obs, "sigma_", sigma, "dist_noise_", dist_noise, "i_",i, ".csv"))[,-1]
# SVCMbeta = t(as.matrix(SVCMbeta))

blue2red_cols = c("#000080", "#000083", "#000087", "#00008B", "#00008F", "#000093", "#000097", "#00009B",
                  "#00009F", "#0000A3", "#0000A7", "#0000AB", "#0000AF", "#0000B3", "#0000B7", "#0000BB",
                  "#0000BF", "#0000C3", "#0000C7", "#0000CB", "#0000CF", "#0000D3", "#0000D7", "#0000DB",
                  "#0000DF", "#0000E3", "#0000E7", "#0000EB", "#0000EF", "#0000F3", "#0000F7", "#0000FB",
                  "#0004FF", "#0008FF", "#000CFF", "#0010FF", "#0014FF", "#0018FF", "#001CFF", "#0020FF",
                  "#0024FF", "#0028FF", "#002CFF", "#0030FF", "#0034FF", "#0038FF", "#003CFF", "#0040FF",
                  "#0044FF", "#0048FF", "#004CFF", "#0050FF", "#0054FF", "#0058FF", "#005CFF", "#0060FF",
                  "#0064FF", "#0068FF", "#006CFF", "#0070FF", "#0074FF", "#0078FF", "#007CFF", "#0080FF",
                  "#0083FF", "#0087FF", "#008BFF", "#008FFF", "#0093FF", "#0097FF", "#009BFF", "#009FFF",
                  "#00A3FF", "#00A7FF", "#00ABFF", "#00AFFF", "#00B3FF", "#00B7FF", "#00BBFF", "#00BFFF",
                  "#00C3FF", "#00C7FF", "#00CBFF", "#00CFFF", "#00D3FF", "#00D7FF", "#00DBFF", "#00DFFF",
                  "#00E3FF", "#00E7FF", "#00EBFF", "#00EFFF", "#00F3FF", "#00F7FF", "#00FBFF", "#00FFFF",
                  "#04FFFB", "#08FFF7", "#0CFFF3", "#10FFEF", "#14FFEB", "#18FFE7", "#1CFFE3", "#20FFDF",
                  "#24FFDB", "#28FFD7", "#2CFFD3", "#30FFCF", "#34FFCB", "#38FFC7", "#3CFFC3", "#40FFBF",
                  "#44FFBB", "#48FFB7", "#4CFFB3", "#50FFAF", "#54FFAB", "#58FFA7", "#5CFFA3", "#60FF9F",
                  "#64FF9B", "#68FF97", "#6CFF93", "#70FF8F", "#74FF8B", "#78FF87", "#7CFF83", "#80FF80",
                  "#83FF7C", "#87FF78", "#8BFF74", "#8FFF70", "#93FF6C", "#97FF68", "#9BFF64", "#9FFF60",
                  "#A3FF5C", "#A7FF58", "#ABFF54", "#AFFF50", "#B3FF4C", "#B7FF48", "#BBFF44", "#BFFF40",
                  "#C3FF3C", "#C7FF38", "#CBFF34", "#CFFF30", "#D3FF2C", "#D7FF28", "#DBFF24", "#DFFF20",
                  "#E3FF1C", "#E7FF18", "#EBFF14", "#EFFF10", "#F3FF0C", "#F7FF08", "#FBFF04", "#FFFF00",
                  "#FFFB00", "#FFF700", "#FFF300", "#FFEF00", "#FFEB00", "#FFE700", "#FFE300", "#FFDF00",
                  "#FFDB00", "#FFD700", "#FFD300", "#FFCF00", "#FFCB00", "#FFC700", "#FFC300", "#FFBF00",
                  "#FFBB00", "#FFB700", "#FFB300", "#FFAF00", "#FFAB00", "#FFA700", "#FFA300", "#FF9F00",
                  "#FF9B00", "#FF9700", "#FF9300", "#FF8F00", "#FF8B00", "#FF8700", "#FF8300", "#FF8000",
                  "#FF7C00", "#FF7800", "#FF7400", "#FF7000", "#FF6C00", "#FF6800", "#FF6400", "#FF6000",
                  "#FF5C00", "#FF5800", "#FF5400", "#FF5000", "#FF4C00", "#FF4800", "#FF4400", "#FF4000",
                  "#FF3C00", "#FF3800", "#FF3400", "#FF3000", "#FF2C00", "#FF2800", "#FF2400", "#FF2000",
                  "#FF1C00", "#FF1800", "#FF1400", "#FF1000", "#FF0C00", "#FF0800", "#FF0400", "#FF0000",
                  "#FB0000", "#F70000", "#F30000", "#EF0000", "#EB0000", "#E70000", "#E30000", "#DF0000",
                  "#DB0000", "#D70000", "#D30000", "#CF0000", "#CB0000", "#C70000", "#C30000", "#BF0000",
                  "#BB0000", "#B70000", "#B30000", "#AF0000", "#AB0000", "#A70000", "#A30000", "#9F0000",
                  "#9B0000", "#970000", "#930000", "#8F0000", "#8B0000", "#870000", "#830000", "#800000")


multi.figs.levelplot = function(zmat, x, y, titles = NULL, cols = blue2red_cols, layout = NULL, xlab = NULL, ylab = NULL, panel = panel.levelplot,
                                colorkey = TRUE) {
  if (is.null(titles)) {
    titles = paste("fig", 1:ncol(zmat))
  }
  pdat = list(z = c(zmat), group = factor(rep(titles, each = nrow(zmat)), levels = titles), x = rep(x, times = ncol(zmat)), y = rep(y, times = ncol(zmat)))
  levelplot(z ~ x + y | group, data = pdat, col.regions = cols, cuts = length(cols) - 1, layout = layout, xlab = xlab,
            ylab = ylab, panel = panel, colorkey = colorkey)
}

show.multi.images.axial = function(imgs,all_coords,zlist,titles,cols=blue2red_cols,col_lim = range(imgs),
                                   layout=NULL,bgcol="black",strip_bgcol = "black",strip_text_col="white",main=""){
  all_coords = as.data.frame(all_coords)
  names(all_coords) = c("x","y","z")
  idx =which(is.element(all_coords$z,zlist))
  
  
  imgs = as.matrix(imgs)
  
  
  
  z = factor(rep(paste("z = ",all_coords$z[idx]),times =ncol(imgs)),levels=paste("z = ",unique(all_coords$z[idx])))
  group = factor(rep(titles,each=length(idx)),levels=titles)
  x = rep(all_coords$x[idx],times=ncol(imgs))
  y = rep(all_coords$y[idx],times=ncol(imgs))
  
  fig=levelplot(imgs[idx,]~x+y | z * group, at=seq(col_lim[1],col_lim[2],length=length(cols)+1),
                col.regions = cols,cuts = length(cols),
                par.settings=list(panel.background=list(col =bgcol)),xlab="x",ylab="y",
                layout=layout,aspect=1.0,main=main,xlim=c(-90, 90),ylim=c(-126,90),
                strip = strip.custom(bg=strip_bgcol,
                                     par.strip.text=list(col=strip_text_col, cex=.8, font=3)))
  return(fig)
}

show.multi.images.axial0 = function(imgs,all_coords,zlist,titles,cols=blue2red_cols,col_lim = range(imgs),
                                    layout=NULL,bgcol="black",strip_bgcol = "black",strip_text_col="white",main=""){
  all_coords = as.data.frame(all_coords)
  names(all_coords) = c("x","y","z")
  idx =which(is.element(all_coords$z,zlist))
  
  
  imgs = as.matrix(imgs)
  
  
  
  z = factor(rep(paste("z = ",all_coords$z[idx]),times =ncol(imgs)),levels=paste("z = ",unique(all_coords$z[idx])))
  group = factor(rep(titles,each=length(idx)),levels=titles)
  x = rep(all_coords$x[idx],times=ncol(imgs))
  y = rep(all_coords$y[idx],times=ncol(imgs))
  
  fig=levelplot(imgs[idx,]~x+y | z * group, at=seq(col_lim[1],col_lim[2],length=length(cols)+1),
                col.regions = cols,cuts = length(cols),
                par.settings=list(panel.background=list(col =bgcol)),xlab="x",ylab="y",
                layout=layout,aspect=1.0,main=main,xlim=c(0, 32),ylim=c(0,32),
                strip = strip.custom(bg=strip_bgcol,
                                     par.strip.text=list(col=strip_text_col, cex=.8, font=3)))
  return(fig)
}


pdf("/Users/yutingduan/Dropbox (University of Michigan)/IOSME/code/different_methods.pdf", height = 20, width = 15)

beta_plot = list()
nr = 1
nc = 3
for(j in 1:3){
  beta_plot[[j]] = show.multi.images.axial0(cbind(true_beta[,j],IOSMEbeta[,j],MUAbeta[,j], SPMbeta[,j],
                                                  BSTbeta[,j],STOREbeta[,j],SVCMbeta[,j], IRRNNbeta[,j]),
                                            voxels,zlist=c(2,4,6,8),
                                            titles = c("True","IOSME", "MUA", "SPM", "BST", "STORE", "SVCM", "IRRNN"),
                                            col_lim = c(-2,2),main = paste("beta",j,sep="_"))
}
plot = arrangeGrob(grobs=beta_plot,nrow=nr,ncol=nc)

grid.arrange(plot, nrow=1,ncol=1)

dev.off()

