# Graphics evaluate densities and coordinates in R; no Python bridge is used.
write_plot <- function(stem,width,height,draw) {
  grDevices::pdf(paste0(stem,".pdf"),width=width,height=height,useDingbats=FALSE)
  tryCatch(draw(),finally=grDevices::dev.off())
  grDevices::png(paste0(stem,".png"),width=width,height=height,units="in",res=180,type="cairo")
  tryCatch(draw(),finally=grDevices::dev.off())
}
make_triangle_figure <- function(g,out="figures") {
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  tr<-g$triangle; means<-tr$means; modes<-tr$mode_search$modes
  xx<-seq(-1.45,4.25,length.out=520); yy<-seq(-1.4,3.9,length.out=490)
  grid<-as.matrix(expand.grid(xx,yy))
  zz<-matrix(exp(mixture_log_density(grid,means)),nrow=length(xx))
  draw<-function() {
    par(mar=c(4.1,4.1,1.1,.8),mgp=c(2.6,.65,0),tcl=-.25,cex=.94)
    contour(xx,yy,zz,levels=seq(.004,.058,length.out=10),col="#3b739b",lwd=1.05,
      xlab="Whitened coordinate 1",ylab="Whitened coordinate 2",asp=1,labcex=.67,bty="l")
    tt<-seq(0,2*pi,length.out=500)
    for(j in 1:3) lines(means[j,1]+.68*cos(tt),means[j,2]+.68*sin(tt),lty=2,col="#67717a",lwd=1.2)
    points(means,pch=4,col="#bc570e",cex=1.15,lwd=1.8)
    points(modes,pch=16,col="#173956",cex=.8)
    legend("topright",c("Component means","Modes","Certified regions"),pch=c(4,16,NA),
      lty=c(NA,NA,2),col=c("#bc570e","#173956","#67717a"),bty="n",cex=.76)
  }
  write_plot(file.path(out,"triangle_modes"),6.4,5.3,draw)
  invisible(list(x=xx,y=yy,density=zz))
}
make_radius_figure <- function(g,out="figures") {
  tab<-g$local_radius$table
  draw<-function() {
    par(mar=c(4,4,1,1),mgp=c(2.5,.7,0))
    b<-barplot(rbind(tab$local_variance,tab$common_variance),beside=TRUE,
      names.arg=paste("Component",1:3),col=c("#315f7c","#d6e2e9"),border=NA,
      ylim=c(0,1.18),ylab="Maximum posterior variance",xlab="Component region")
    abline(h=.9,lty=2,col="#bc570e"); abline(h=1,lty=3,col="grey45")
    legend("topleft",c("Separate radii","Common radius",expression(beta==0.9)),
      fill=c("#315f7c","#d6e2e9",NA),border=NA,lty=c(NA,NA,2),col=c(NA,NA,"#bc570e"),bty="o",bg="white",box.col="white",cex=.83)
  }
  write_plot(file.path(out,"local_radii"),6.4,4.4,draw)
}
make_fit_figure <- function(tab,initial_ll,out="figures") {
  fixed<-tab[tab$check=="fixed",]; joint<-tab[tab$check=="joint",]
  draw<-function() {
    par(mfrow=c(1,2),mar=c(4.1,4.3,1.5,.7),mgp=c(2.8,.7,0),cex=.86)
    plot(c(0,fixed$step),c(initial_ll,fixed$log_likelihood),type="b",pch=16,
      col="#315f7c",xlab="Accepted fixed-budget step",ylab="Mean log likelihood",xaxt="n",bty="l")
    axis(1,at=0:3); title("Fixed-budget trajectory",cex.main=1)
    vals<-rbind(tab$likelihood_gain,tab$auxiliary_gain,tab$proximal_penalty)
    barplot(vals,beside=TRUE,names.arg=c("Fixed 1","Fixed 2","Fixed 3","Joint 1"),
      col=c("#315f7c","#90afc2","#d6e2e9"),border=NA,ylim=c(0,max(vals)*1.15),ylab="Per-observation gain or penalty",bty="l")
    legend("topleft",c("Log likelihood","EM auxiliary","Proximal penalty"),
      fill=c("#315f7c","#90afc2","#d6e2e9"),border=NA,bty="n",cex=.72)
  }
  write_plot(file.path(out,"fitting_checks"),9,4.1,draw)
}
