load("summariesv9.RData")
models <- names(sumup)
models.name <- c("1 main", "2 main", "3 main", "1 main X", "2 main X", "2 main 1 int", "2 main 1 int AX")
postscript("../Figs/fig_sim.eps", height=8, width=7, pointsize=12, onefile=FALSE, paper="special", horizontal=FALSE)
#par(mfrow=c(3,2), las=1, mar=c(7.1,5.1,2.6,3.1))
par(mfrow=c(2,2), las=1)
par(oma = c(2, 0, 0, 0))
par(mar = c(6.5, 5, 2, 2)) # make the plots be closer together
#Fig1
typeI <- vector()
for (model.id in models){
  typeI[[model.id]] <- 1000 - attr(sumup[[model.id]][[1]],"qtl.extra")[1]
  }
x <- 1:length(typeI)-0.5
plot(0,0,type="n", main="A", xlab="", ylab="Main effect Type I error (%)", ylim=c(0,25),
     xlim=c(0.5, length(typeI)+0.5), xaxt="n", yaxs="i")
abline(h=seq(10,20,by=5), lty=2, col="gray")
abline(h=5, lty=2, col="red")
#methods <- c("1.42.8","3MBest","H-L")
blue <- "slateblue"
red <- "violetred"
black <- "black"
green <- broman::brocolors("web")["green"]
thecol <- c(black, green)
for(j in 1:2) {
  typeI <- vector()
  for (model.id in models){
    typeI[[model.id]] <- 1000 - attr(sumup[[model.id]][[j]],"qtl.extra")[1]
    }
  x <- 1:length(typeI)-(2-j)*0.2

  segments(x-0.1, typeI/10,
           x+0.1, typeI/10, lwd=2,col=thecol[j]
           )
  ci <- matrix(ncol=2, nrow=length(typeI))
  for(k in 1:length(typeI))
    ci[k,] <- binom.test(typeI[k], 1000)$conf.int*100
  segments(x, ci[,1], x, ci[,2],col=thecol[j])

  axis(side=1, at=1:length(typeI), tick=FALSE, lab=rep("", length(typeI)))
  for(k in 1:length(typeI))
    axis(side=1, at=k, models.name[k], las=2,
          tick=FALSE)
#  mtext(side=3, line=1, methods[j])
  #mtext(side=1, line=6, "Generating model", cex=1.2)
}
legend("topleft",legend=c("XeqA","XneA"),fill=thecol, bg="white")

#Fig2
library(qtl)
data(map10)
LA <- sum(chrlen(map10[-20]))
LX <- sum(chrlen(map10[20]))
n.sim <- 1000
x <- 1:length(models)-0.5
plot(0,0,type="n", main="B", xlab="",
     ylab=expression(paste(log[2], " ratio of Type I error in Autosomes vs X Chr")),
     ylim=c(-4,4),
     xlim=c(0.5, length(models)+0.5), xaxt="n", yaxs="i")
abline(h=seq(-3,3,by=1), lty=2, col="gray")
abline(h=0, lty=2, col="red")

for(j in 1:2) {
  typeI.A <- vector()
  typeI.X <- vector()
  for (model.id in models){
    typeI.A[[model.id]] <- 1000 - attr(sumup[[model.id]][[j]],"qtl.extra.A")[1]
    typeI.X[[model.id]] <- 1000 - attr(sumup[[model.id]][[j]],"qtl.extra.X")[1]
  }
  typeI <- log2((typeI.A/LA)/(typeI.X/LX))
  print(typeI)
  x <- 1:length(typeI)-(2-j)*0.2

  segments(x-0.1, typeI,
           x+0.1, typeI, lwd=2,col=thecol[j]
  )
  ci <- matrix(ncol=2, nrow=length(typeI))
  typeI.A.p <- typeI.A/n.sim
  typeI.X.p <- typeI.X/n.sim
  for(k in 1:length(typeI)) {
    # With A+X=alpha assumption
    sd1 <- sqrt(abs((1/n.sim)*((1-typeI.X.p[k])/typeI.X.p[k]-1*(1-typeI.A.p[k])/typeI.A.p[k])))
    # use upper bound of the standard error derived from Cauchy Schwartz inequality: Var(X-Y)<=(sigma_x+sigma_y)^2
    sd2 <- sqrt((1/n.sim)*((1-typeI.X[k]/n.sim)/(typeI.X[k]/n.sim)))+sqrt((1/n.sim)*((1-typeI.A[k]/n.sim)/(typeI.A[k]/n.sim)))
    # with independent assumption
    sd3 <- sqrt((1/n.sim)*((1-typeI.X[k]/n.sim)/(typeI.X[k]/n.sim))+(1/n.sim)*((1-typeI.A[k]/n.sim)/(typeI.A[k]/n.sim)))
    # using emperical covariance of 2 Bernoulli distribution and Taylor approximation of var(f(x,y)) with f(x,y)=log(x/y)
    sd4 <- sqrt(((1-typeI.X.p[k])/(n.sim*typeI.X.p[k]))+((1-typeI.A.p[k])/(n.sim*typeI.A.p[k]))-2*n.sim*(cov(attr(sumup[[k]][[j]],"qtl.extra.A.vec"),attr(sumup[[k]][[j]],"qtl.extra.X.vec")))/(typeI.X[k]*typeI.A[k]))

    #     print(sd)
    sd <- sd4
    ci[k,1] <- typeI[k]-2*sd
    ci[k,2] <- typeI[k]+2*sd
  }
  segments(x, ci[,1], x, ci[,2],col=thecol[j])

  axis(side=1, at=1:length(typeI), tick=FALSE, lab=rep("", length(typeI)))
  for(k in 1:length(typeI))
    axis(side=1, at=k, models.name[k], las=2,
         tick=FALSE)
  #mtext(side=1, line=6, "Generating model", cex=1.2)
}
legend("topright", legend=c("XeqA","XneA"),fill=thecol, bg="white")

#Fig3
typeI <- vector()
for (model.id in models){
  typeI[[model.id]] <- 1000 - attr(sumup[[model.id]][[j]],"int.extra")[1]
}
x <- 1:length(typeI)-0.5
plot(0,0,type="n", main="C", xlab="", ylab="Interaction Type I error (%)", ylim=c(0,17),
     xlim=c(0.5, length(typeI)+0.5), xaxt="n", yaxs="i")
abline(h=seq(2,16,by=2), lty=2, col="gray")
abline(h=5, lty=2, col="red")

for(j in 1:2) {
  typeI <- vector()
  for (model.id in models){
    typeI[[model.id]] <- 1000 - attr(sumup[[model.id]][[j]],"int.extra")[1]
  }
  x <- 1:length(typeI)-(2-j)*0.2

  segments(x-0.1, typeI/10,
           x+0.1, typeI/10, lwd=2,col=thecol[j]
  )
  ci <- matrix(ncol=2, nrow=length(typeI))
  for(k in 1:length(typeI))
    ci[k,] <- binom.test(typeI[k], 1000)$conf.int*100
  segments(x, ci[,1], x, ci[,2],col=thecol[j])

  axis(side=1, at=1:length(typeI), tick=FALSE, lab=rep("", length(typeI)))
  for(k in 1:length(typeI))
    axis(side=1, at=k, models.name[k], las=2,
         tick=FALSE)
  #mtext(side=1, line=6, "Generating model", cex=1.2)
}
legend("topleft",legend=c("XeqA","XneA"),fill=thecol, bg="white")

#Fig4
plot.value <- vector()
plot.feature <- "qtl.missed"
for (model.id in models){
  plot.value[[model.id]] <- attr(sumup[[model.id]][[1]],plot.feature)[1]
}
x <- 1:length(plot.value)-0.5
plot(0,0,type="n", main="D", xlab="", ylab="Power to detect all QTL (%)", ylim=c(0,100),
     xlim=c(0.5, length(plot.value)+0.5), xaxt="n", yaxs="i")
abline(h=seq(10,90,by=10), lty=2, col="gray")
# abline(h=5, lty=2, col="red")

for(j in 1:2) {
  plot.value <- vector()
  for (model.id in models){
    plot.value[[model.id]] <-attr(sumup[[model.id]][[j]],plot.feature)[1]
  }
  x <- 1:length(plot.value)-(2-j)*0.2

  segments(x-0.1, plot.value/10,
           x+0.1, plot.value/10, lwd=2,col=thecol[j]
  )
  ci <- matrix(ncol=2, nrow=length(plot.value))
  for(k in 1:length(plot.value))
    ci[k,] <- binom.test(plot.value[k], 1000)$conf.int*100
  segments(x, ci[,1], x, ci[,2],col=thecol[j])

  axis(side=1, at=1:length(plot.value), tick=FALSE, lab=rep("", length(plot.value)))
  for(k in 1:length(plot.value))
    axis(side=1, at=k, models.name[k], las=2,
         tick=FALSE)
  #mtext(side=1, line=6, "Generating model", cex=1.2)
}
legend("bottomleft",legend=c("XeqA", "XneA"),fill=thecol, bg="white")

#mtext('Generating model', side = 1, outer = TRUE, line = 1)
dev.off()
