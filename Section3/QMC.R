##############  this file is to plot Figgire 1 to illustrate the process of RQMC

library(randtoolbox)

QMC <- sobol(200, dim = 2, scrambling = 3)

uu <- c(0.107036, -0.2012338)

QMC1 <- sweep(QMC, 2, uu, "+")

QMC <- as.data.frame(QMC)
names(QMC) <- c("x","y")

QMC1 <- as.data.frame(QMC1)
names(QMC1) <- c("x","y")

QMC1_in <- QMC1[which(QMC1$x<1&QMC1$y>0),]
QMC1_out <- QMC1[-which(QMC1$x<1&QMC1$y>0),]

randomization <- function(QMC, lower = rep(0, ncol(QMC)), upper = rep(1, ncol(QMC)))
{
  n <- nrow(QMC)
  p <- ncol(QMC)
  lowerm <- matrix(rep(lower,rep(n,p)),n,p)
  upperm <- matrix(rep(upper,rep(n,p)),n,p)
  right <- QMC - upperm
  QMC[right>0] <- lowerm[right>0] + right[right>0]
  left <- lowerm - QMC
  QMC[left>0] <- upperm[left>0] - left[left>0]

  return(QMC)
}

lower <- rep(0,2)
upper <- rep(1,2)

RQMC <- randomization(QMC1, lower = lower, upper = upper)

RQMC <- as.data.frame(RQMC)
names(RQMC) <- c("x","y")

RQMC1 <- RQMC[-which(QMC1$x<1&QMC1$y>0),]


library(ggplot2)

ggplot(data = QMC, aes(x = x, y = y)) + geom_point() +
  xlim(0, 1.15) + ylim(-0.25, 1) + 
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha =.1, fill = "#006699") + theme_bw()+
  theme(axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank())

ggplot() + geom_point(data = QMC1_out, aes(x = x, y = y), shape = 17, col = "red") +
  geom_point(data = QMC1_in, aes(x = x, y = y)) +
  xlim(0, 1.15) + ylim(-0.25, 1) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha =.1, fill = "#006699") + theme_bw() +
  theme(axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank())



ggplot() + geom_point(data = RQMC1, aes(x = x, y = y), shape = 17, col = "red") +
  geom_point(data = QMC1_in, aes(x = x, y = y)) +
  xlim(0, 1.15) + ylim(-0.25, 1) +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha =.1, fill = "#006699") + theme_bw()+
  theme(axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank())



  



