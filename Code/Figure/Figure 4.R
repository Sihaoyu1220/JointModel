chippy4 <- function(time, T=80, P1=11, lambda1 = 0.2, P2 = 2.5, beta1 = 1.8, beta2=8, beta3 = 0.3, beta4 = 1.1) log10( exp(P1-lambda1*time) + exp(P2))*(time<T)+ (beta1 * (time-T) / ((time-T)  + exp(beta2-beta3*(time-T))) + beta4 )*(time>=T)
curve(chippy4,1, 140, ylab = expression(paste("Viral load (", log[10], " -Transformed)")) , xlab="Day", xaxt="n",cex.lab=1.3, mgp=c(2,0.5,0))
axis(1, at = seq(10, 130, by = 20))
