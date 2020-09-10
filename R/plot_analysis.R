library(qtl)
load("analysis.RData")
#Load back the probe stepwise result
library(knitr)
kable(list(summary(out_10003837305), summary(out_10003837305_curr)), format = "latex")
green <- broman::brocolors("web")["green"]
black <- "black"
postscript("../Figs/fig_prob305.eps", height=5, width=7, pointsize=12, onefile=FALSE, paper="special", horizontal=FALSE)
plotLodProfile(out_10003837305_curr, col=black, qtl.labels=FALSE, ylab="Profile LOD score")
abline(h=seq(10, 40, by=10), lty=2, col="gray")
plotLodProfile(out_10003837305_curr, col=black, qtl.labels=FALSE, add=TRUE)
plotLodProfile(out_10003837305, ,add=TRUE, col=green, qtl.labels=FALSE)
legend("topleft",legend=c("XeqA", "XneA"),fill = c(black, green), bg="white")
dev.off()
