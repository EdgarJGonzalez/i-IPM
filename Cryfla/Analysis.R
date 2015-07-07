# FITTING AN IPM TO Roberto Salguero-Gomez DATA

require(lme4)
require(MASS)
require(scales)

get.binCI <- function(X,N) {
 CI <- matrix(ncol = 2, nrow = length(X))
 for (i in 1:length(X)) {
  if(N[i] == 0)
   CI[i, ] = rbind(0, 0)
  if(N[i] != 0)
   CI[i, ] = rbind(setNames(c(binom.test(X[i],N[i])$conf.int),c("lwr","upr")))
 }
 return(CI)
}

setwd("/Users/Edgar/Documents/UNAM/Postdoc/Ben/Paper1/Cryfla/Data")

data <- read.csv("Data.txt", sep = "\t")
# removing indivuals subject to experimental drought, individuals with size 0
#  and year 2002 because no data were recorded on that year
data <- subset(data, Treatment == "C" & Year != 2002)
data <- data[-which(data$Size == 0), ]
# scaling Size
mu.Size <- mean(log10(data$Size), na.rm = T)
sd.Size <- sd(log10(data$Size), na.rm = T)
data <- transform(data, Plot = as.factor(Plot), Quadrat = as.factor(Quadrat), 
 log.Size = ifelse(is.na(Size), NA, scale(log10(Size), scale = sd.Size)),
 sqrt.Size = ifelse(is.na(Size), NA, scale(sqrt(Size), scale = sd.Size)))
palette <- data.frame(Year = unique(data$Year), Year.col =
 sort(rainbow(length(unique(data$Year)))), Year.lty = rep(c(1, 2, 3), 
 length.out = length(unique(data$Year))))
data <- merge(data, palette)

### DATA

## Survival data
data.s <- c()
for (ind in levels(data$ID)) {
 data.ind <- subset(data, data$ID == ind & is.na(data$log.Size) != T)
 if(length(data.ind$Year) == 0)
  next
 data.ind <- transform(data.ind, Survival = 0)
 out <- c()
 if (length(data.ind$Year) == 1) {
  if (data.ind$Year == 2001 | data.ind$Year == 2012)
   next
  if (data.ind$Year != 2001 & data.ind$Year != 2012) {
   data.ind$Survival <- 0
   data.ind$Year <- data.ind$Year + 1
   data.s <- rbind(data.s, data.ind) 
   next
  }
 }
 if (length(data.ind$Year) > 1) {
  if (any(data.ind$Year == 2001))
   out <- c(out, which(data.ind$Year == 2001))
  if (any(data.ind$Year == 2012))
   out <- c(out, which(data.ind$Year == 2012))
  if(length(out) > 0)
   data.ind <- data.ind[-out, ]
 }
 if (length(data.ind$Year) == 0)
  next
 if (length(data.ind$Year) == 1) {
  data.ind$Survival <- 0
  data.ind$Year <- data.ind$Year + 1
  data.s <- rbind(data.s, data.ind)    
  next
 }
 if (length(data.ind$Year) > 1) {
  for (i in 1:(length(data.ind$Year)-1)) {
   data.ind$Survival[i] <- 1
   data.ind$Year[i] <- data.ind$Year[i] + 1
  }
  data.ind$Survival[length(data.ind$Year)] <- 0
  data.ind$Year[length(data.ind$Year)] <- data.ind$Year[length(data.ind$Year)] + 1
  data.s <- rbind(data.s, data.ind)  
 }
} 
data.s <- droplevels(data.s)

## Growth data
data.G <- subset(data, is.na(data$log.Size) == F)
data.g <- c()
for (ind in levels(data.G$ID)) {
 data.ind <- subset(data.G, data.G$ID == ind)
 if (dim(data.ind)[1] > 1) {
  log.size1 <- data.ind$log.Size[1:(length(data.ind$log.Size)-1)]
  log.size2 <- data.ind$log.Size[2:length(data.ind$log.Size)]
  sqrt.size1 <- data.ind$sqrt.Size[1:(length(data.ind$sqrt.Size)-1)]
  sqrt.size2 <- data.ind$sqrt.Size[2:length(data.ind$sqrt.Size)]
  data.ind <- data.ind[-1,]
  data.ind <- transform(data.ind, log.Size1 = log.size1, 
   log.Size2 = log.size2, sqrt.Size1 = sqrt.size1, 
   sqrt.Size2 = sqrt.size2)
  data.g <- rbind(data.g, data.ind)
 }
}
data.g <- subset(data.g, data.g$Year != 2003)

## Fecundity data
data.f <- subset(data, is.na(data$Fert) == F)

## Seedling data
data.ns <- c()
for (ind in levels(data$ID)) {
 data.ind <- subset(data, data$ID == ind)
 if (dim(data.ind)[1] > 1)
  for (i in 2:length(unique(data.ind$Year)))
   if (is.na(data.ind$log.Size[i-1]) & !is.na(data.ind$log.Size[i]))
 data.ns <- rbind(data.ns, data.ind[i, ])
}
# individuals that are already reproducing on the first year they're recorded 
#  are not considered seedlings
data.ns <- subset(data.ns, data.ns$Fert == 0)

## Data visualization
x <- seq(min(data$Size, na.rm = T), max(data$Size, na.rm = T), 1)
hist(data$Size, breaks = x, main = "Size distribution", xlab = "Size")
hist(data.f$Size[which(data.f$Fert > 0)], breaks = x, add = T, border = "red")
hist(data.ns$Size, breaks = x, add = T, border = "blue")
legend.txt <- c("all", "reproductive", "newborns")
legend("topright",legend.txt, text.col = c("black", "red", "blue"), bty = "n")

inc <- (max(data$log.Size, na.rm = T)-min(data$log.Size, na.rm = T))/50
x <- seq(min(data$log.Size, na.rm = T), max(data$log.Size, na.rm = T), inc)
hist(data$log.Size, breaks = x, main = "log10(Size) distribution", 
 xlab = "scale(log10(Size))")
hist(data.f$log.Size[which(data.f$Fert > 0)], breaks = x, add = T, 
 border = "red")
hist(data.ns$log.Size, breaks = x, add = T, border = "blue")
legend.txt <- c("all", "reproductive", "newborns")
legend("topright",legend.txt, text.col = c("black", "red", "blue"), bty = "n")

inc <- (max(data$sqrt.Size, na.rm = T)-min(data$sqrt.Size, na.rm = T))/50
x <- seq(min(data$sqrt.Size, na.rm = T), max(data$sqrt.Size, na.rm = T), inc)
hist(data$sqrt.Size, breaks = x, main = "sqrt(Size) distribution", 
 xlab = "scale(sqrt(Size))")
hist(data.f$sqrt.Size[which(data.f$Fert > 0)], breaks = x, add = T, 
 border = "red")
hist(data.ns$sqrt.Size, breaks = x, add = T, border = "blue")
legend.txt <- c("all", "reproductive", "newborns")
legend("topright",legend.txt, text.col = c("black", "red", "blue"), bty = "n")


### ANALYSES

## Survival analysis

modS0 <- glm(Survival ~ poly(log.Size,2), family = binomial, data = data.s)
# AIC = 4114
# (Intercept)  poly(log.Size, 2)1  poly(log.Size, 2)2  
# 0.5954       25.5139             -6.5696  

modS1 <- glmer(Survival ~ poly(log.Size,2) + (1|Year), family = binomial, 
 data = data.s)
# AIC = 3497.351
# Fixed effects:
#  (Intercept)  poly(log.Size, 2)1  poly(log.Size, 2)2  
#  0.4595       32.6414             -7.2042 
# Random effects:
#  (Intercept)
#  1998  0.93444658
#  1999  0.07147408
#  2000  1.23614464
#  2001 -1.53157130
#  2004  0.72446256
#  2005  1.07848759
#  2006  0.40153117
#  2007 -0.14398483
#  2008  0.20967902
#  2009  0.58294929
#  2010 -0.19040680
#  2011  0.81136448
#  2012 -4.13401881

# plots modS0
{
 # fitting
 
 data.s0 <- subset(data.s, Survival == 0)
 counts.0 <- table(data.s0$log.Size)
 plot(as.numeric(names(counts.0)), rep(0, length(counts.0)), 
  bty = "n", xlim = c(min(x), max(x)), ylim = c(0, 1), yaxs = "i", 
  main = "Survival", ylab = "Survival probability", 
  xlab = "scale(log10(Size1))", col = "grey", pch = 16, 
  cex =  sqrt(counts.0/pi))
 data.s1 <- subset(data.s, Survival == 1)
 counts.1 <- table(data.s1$log.Size)
 points(as.numeric(names(counts.1)), rep(1, length(counts.1)), col = "grey",
  pch = 16, cex = sqrt(counts.1/pi))
 lines(x, predict(modS0, newdata = list(log.Size = x), type = "response"), 
  lwd = 4)
 get.binCI(S0.counts,S0.counts+S1.counts)
 
 # residuals
 residuals <- residuals(modS0, "deviance")
 data.s <- transform(data.s, sum.res.0 = log.Size+residuals)
 counts <- table(data.s$sum.res.0)
 counts.col <- rep(1, dim(data.s)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.s)[1])
   if (round(data.s$sum.res.0[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modS0, type = "response")
 plot(fitted, residuals, ylim = c(-4, 4), pch = 16, bty = "n", col = "grey", 
  cex =  sqrt(counts.col/pi), main = "Survival model residuals", 
  xlab = "scale(log10(Size))", ylab = "Deviance residuals")
 abline(0, 0)
}

# plots modS1
{
 # fitting
 plot(as.numeric(names(counts.0)), rep(0, length(counts.0)), 
  bty = "n", xlim = c(min(x), max(x)), ylim = c(0, 1), yaxs = "i", 
  main = "Survival", ylab = "Survival probability", 
  xlab = "scale(log10(Size1))", col = "grey", pch = 16, 
  cex =  sqrt(counts.0/pi))
 points(as.numeric(names(counts.1)), rep(1, length(counts.1)), col = "grey",
  pch = 16, cex = sqrt(counts.1/pi))
 lines(x, predict(modS1, newdata = list(log.Size = x), type = "response", 
  re.form = ~ 0), lwd = 4)
 coef <- coef(modS1)$Year
 rownames(coef) <- as.numeric(rownames(coef(modS1)$Year))
 poly <- poly(data.s$log.Size,2)
 legend <- data.frame(Year = unique(data.s$Year), 
  Year.col = rep(NA,length(unique(data.s$Year))), 
  Year.lty = rep(NA,length(unique(data.s$Year))))
 class(legend$Year.col) <- "character"
 for (i in 1:length(unique(data.s$Year))) {
  year.i <- unique(data.s$Year)[i]
  subset <- subset(data.s, data.s$Year == year.i)
  subset <- droplevels(subset)
  coef.year <- coef[which(names(coef) == year.i)]
  x <- seq(min(subset$log.Size), max(subset$log.Size), 0.1)
  X <- predict(poly, x)
  coefs <- coef[which(rownames(coef) == year.i), ]
  y <- plogis(as.numeric(coefs[1]) + as.numeric(coefs[2])*X[, 1] 
   + as.numeric(coefs[3])*X[, 2])
  legend$Year.col[i] <- as.character(subset$Year.col[1])
  legend$Year.lty[i] <- subset$Year.lty[1]
  lines(x, y, col = alpha(legend$Year.col[i], 0.5), 
   lty = legend$Year.lty[i], lwd = 2)
  #lines(x[-length(x)], y[-length(y)], col = subset$Year.col[1])
  # arrows(x[length(x)-1], y[length(y)-1], x[length(x)], y[length(y)], 
  # col = subset$Year.col[1], length = 0.1)
 }
 legend <- legend[order(legend$Year),]
 legend("bottomright", legend = c("All",legend[, 1]), col = c("black", 
  alpha(legend[, 2], 0.5)), lty = c(1, legend[, 3]), bty = "n", 
  lwd = c(4, rep(2, length(legend$Year))))
    
 # residuals
 residuals <- residuals(modS1, "deviance")
 data.s <- transform(data.s, sum.res.1 = log.Size+residuals)
 counts <- table(data.s$sum.res.1)
 counts.col <- rep(1, dim(data.s)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.s)[1])
   if (round(data.s$sum.res.1[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modS1, type = "response")
 plot(fitted, residuals, col = "grey", ylim = c(-4, 4), 
  pch = 16, bty = "n", cex =  sqrt(counts.col/pi), 
  main = "Survival model residuals", xlab = "scale(log10(Size))", 
  ylab = "Deviance residuals")  
 abline(0, 0)
}

## Growth analysis

modG0 <- lm(log.Size2 ~ log.Size1, data = data.g)
# AIC = 4515.866
# (Intercept)    log.Size1  
# 0.1927         0.6445
log10(summary(modG0)$sigma) # -0.1928407 

modG1 <- lmer(log.Size2 ~ log.Size1 + (1|Year), data = data.g, REML = F) 
# AIC = 4309.704
# Fixed effects:
#  (Intercept)    log.Size1  
#  0.2397         0.6410
log10(sigma(modG1)) # -0.2160386
# Random effects:
#  (Intercept)
#  1998 -0.16299298
#  1999 -0.35743481
#  2000  0.10390452
#  2001  0.02953927
#  2004  0.33010814
#  2005  0.10099882
#  2006 -0.28930445
#  2007 -0.23225019
#  2008 -0.03121846
#  2009  0.20433340
#  2010  0.13130298
#  2011  0.14270514
#  2012  0.03030862

# plots modG0
{
 # fitting
 data.g <- transform(data.g, sum.Size = log.Size1+log.Size2)
 counts <- table(data.g$sum.Size)
 counts.col <- rep(1, dim(data.g)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.g)[1])
   if (round(data.g$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 x <- seq(min(min(data.g$log.Size1), min(data.g$log.Size2)), 
  max(max(data.g$log.Size1), max(data.g$log.Size2)), 0.1)
 plot(data.g$log.Size1, data.g$log.Size2, bty = "n", xlim = c(min(x), max(x)), 
  ylim = c(min(x), max(x)), main = "Growth", ylab = "scale(log10(Size2))", 
  xlab = "scale(log10(Size1))", col = "grey", pch = 16, cex =  sqrt(counts.col/pi))
 lines(x, predict(modG0, newdata = list(log.Size1 = x), type = "response"), 
  lwd = 4)
 
 # residuals
 residuals <- residuals(modG0, "deviance")
 data.g <- transform(data.g, sum.res.0 = log.Size1+residuals)
 counts <- table(data.g$sum.res.0)
 counts.col <- rep(1, dim(data.g)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.g)[1])
   if (round(data.g$sum.res.0[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modG0, type = "response")
 plot(fitted, residuals, ylim = c(-4, 4), pch = 16, bty = "n", col = "grey", 
  cex =  sqrt(counts.col/pi), main = "Growth model residuals", 
  xlab = "scale(log10(Size1))", ylab = "Deviance residuals")
 abline(0, 0)
}

# plots modG1
{
 # fitting
 data.g <- transform(data.g, sum.Size = log.Size1+log.Size2)
 counts <- table(data.g$sum.Size)
 counts.col <- rep(1, dim(data.g)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.g)[1])
   if (round(data.g$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 plot(data.g$log.Size1, data.g$log.Size2, bty = "n", xlim = c(min(x), max(x)), 
  ylim = c(min(x), max(x)), main = "Growth", ylab = "scale(log10(Size2))", 
  xlab = "scale(log10(Size1))", col = "grey", pch = 16, cex =  sqrt(counts.col/pi))
 ranef <- ranef(modG1)$Year
 rownames(ranef) <- as.numeric(rownames(ranef(modG1)$Year))
 legend <- data.frame(Year = unique(data.g$Year), 
  Year.col = rep(NA,length(unique(data.g$Year))), 
  Year.lty = rep(NA,length(unique(data.g$Year))))
 class(legend$Year.col) <- "character"
 for (i in 1:length(unique(data.g$Year))) {
  year.i <- unique(data.g$Year)[i]
  subset <- subset(data.g, data.g$Year == year.i)
  subset <- droplevels(subset)
  ranef.year <- ranef[which(as.numeric(rownames(ranef)) == year.i), 1]
  x <- seq(min(min(subset$log.Size1), min(subset$log.Size2)), 
   max(max(subset$log.Size1), max(subset$log.Size2)), 0.1)
  y <- predict(modG1, newdata = list(log.Size1 = x), re.form = ~ 0) + ranef.year
  legend$Year.col[i] <- as.character(subset$Year.col[1])
  legend$Year.lty[i] <- subset$Year.lty[1]
  lines(x, y, col = alpha(legend$Year.col[i], 0.5), 
   lty = legend$Year.lty[i], lwd = 2)
  #lines(x[-length(x)], y[-length(y)], col = subset$Year.col[1])
  # arrows(x[length(x)-1], y[length(y)-1], x[length(x)], y[length(y)], 
  # col = subset$Year.col[1], length = 0.1)
 }
 lines(x, predict(modG1, newdata = list(log.Size1 = x), type = "response", 
  re.form = ~ 0), lwd = 4)  
 legend <- legend[order(legend$Year),]
 legend("bottomright", legend = c("All",legend[, 1]), col = c("black", 
  alpha(legend[, 2], 0.5)), lty = c(1, legend[, 3]), bty = "n", 
  lwd = c(4, rep(2, length(legend$Year))))
    
 # residuals
 residuals <- residuals(modG1, "deviance")
 data.g <- transform(data.g, sum.res.1 = log.Size1+residuals)
 counts <- table(data.g$sum.res.1)
 counts.col <- rep(1, dim(data.g)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.g)[1])
   if (round(data.g$sum.res.1[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modG1, type = "response")
 plot(fitted, residuals, col = "grey", ylim = c(-4, 4), 
  pch = 16, cex =  sqrt(counts.col/pi), bty = "n", main = "Growth model residuals", 
  xlab = "scale(log10(Size1))", ylab = "Deviance residuals")
 abline(0, 0)
}

## Fecundity analysis

# Poisson

modF0 <- glm(Fert ~ log.Size, family = poisson, data = data.f)
# AIC = 11643
# (Intercept)     log.Size  
# -0.328          1.415  
log10(fixef(modF0)[2]) # 0.1508113


modF1 <- glmer(Fert ~ log.Size + (1|Year),  family = poisson, data = data.f)
# AIC = 10369.564
# Fixed effects:
#  (Intercept)   log.Size  
#  -0.6232       1.4185
log10(fixef(modF1)[2]) # 0.1518317 
# Random effects:
#  (Intercept)
#  1997  0.35488996
#  1998  0.73265521
#  1999 -0.29343535
#  2000  0.54886230
#  2001  0.50634282
#  2003 -2.14841165
#  2004  0.28298810
#  2005  0.68394699
#  2006 -0.32082356
#  2007 -0.37521829
#  2008 -0.27206480
#  2009  0.15783656
#  2010  0.07334472
#  2011  0.28812358
#  2012 -0.17826030

# plots modF0
{
 # fitting
 data.f <- transform(data.f, sum.Size = log.Size+Fert)
 counts <- table(data.f$sum.Size)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 x <- seq(min(data.f$log.Size, na.rm = T), max(data.f$log.Size, na.rm = T), 0.1)
 plot(data.f$log.Size, data.f$Fert, pch = 16, col = "grey", 
  cex = sqrt(counts.col/pi), xlim = c(min(x), max(x)), main = "Fecundity", 
  ylab = "Number of flowers", xlab = "log.Size", bty = "l")
 lines(x, predict(modF0, newdata = list(log.Size = x), type = "response"), 
  lwd = 4)
  
 # residuals
 residuals <- residuals(modF0, "deviance")
 data.f <- transform(data.f, sum.res.0 = log.Size+residuals)
 counts <- table(data.f$sum.res.0)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.res.0[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modF0, type = "response")
 plot(fitted, residuals, ylim = c(-4, 4), pch = 16, bty = "n", col = "grey", 
  cex = sqrt(counts.col/pi), main = "Fecundity model residuals", 
  xlab = "scale(log10(Size))", ylab = "Deviance residuals")
 abline(0, 0)
}

# plots modF1
{
 # fitting
 data.f <- transform(data.f, sum.Size = log.Size+Fert)
 counts <- table(data.f$sum.Size)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 plot(data.f$log.Size, data.f$Fert, bty = "n", xlim = c(min(x), max(x)), 
  main = "Fecundity", ylab = "Number of flowers", 
  xlab = "scale(log10(Size))", col = "grey", pch = 16, 
  cex = sqrt(counts.col/pi))
 ranef <- ranef(modF1)$Year
 rownames(ranef) <- as.numeric(rownames(ranef(modF1)$Year))
 legend <- data.frame(Year = unique(data.f$Year), 
  Year.col = rep(NA,length(unique(data.f$Year))), 
  Year.lty = rep(NA,length(unique(data.f$Year))))
 class(legend$Year.col) <- "character"
 for (i in 1:length(unique(data.f$Year))) {
  year.i <- unique(data.f$Year)[i]
  subset <- subset(data.f, data.f$Year == year.i)
  subset <- droplevels(subset)
  ranef.year <- ranef[which(as.numeric(rownames(ranef)) == year.i), 1]
  x <- seq(min(subset$log.Size), max(subset$log.Size), 0.1)
  y <- exp(predict(modF1, newdata = list(log.Size = x), re.form = ~ 0) + 
   ranef.year)
  legend$Year.col[i] <- as.character(subset$Year.col[1])
  legend$Year.lty[i] <- subset$Year.lty[1]
  lines(x, y, col = alpha(legend$Year.col[i], 0.5), 
   lty = legend$Year.lty[i], lwd = 2)
  #lines(x[-length(x)], y[-length(y)], col = subset$Year.col[1])
  # arrows(x[length(x)-1], y[length(y)-1], x[length(x)], y[length(y)], 
  # col = subset$Year.col[1], length = 0.1)
 }
 x <- seq(min(data.f$log.Size, na.rm = T), max(data.f$log.Size, na.rm = T), 
  0.1)
 lines(x, predict(modF1, newdata = list(log.Size = x), type = "response", 
  re.form = ~ 0), lwd = 4)  
 legend <- legend[order(legend$Year),]
 legend("topleft", legend = c("All",legend[, 1]), col = c("black", 
  alpha(legend[, 2], 0.5)), lty = c(1, legend[, 3]), bty = "n", 
  lwd = c(4, rep(2, length(legend$Year))))
    
 # residuals
 residuals <- residuals(modF1, "deviance")
 data.f <- transform(data.f, sum.res.1 = log.Size+residuals)
 counts <- table(data.f$sum.res.1)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.res.1[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modF1, type = "response")
 plot(fitted, residuals, col = "grey", ylim = c(-4, 4), 
  pch = 16, cex =  sqrt(counts.col/pi), bty = "n", 
  main = "Fecundity model residuals", 
  xlab = "scale(log10(Size))", ylab = "Deviance residuals")
 abline(0, 0)
}

# Neg-binomial

modF0 <- glm.nb(Fert ~ log.Size, data = data.f)
# AIC = 9517.5
# (Intercept)   log.Size  
# -0.4915       1.6483 

modF1 <- glmer.nb(Fert ~ log.Size + (1|Year), data = data.f)
# AIC = 11860.38
# Fix effects:
#  (Intercept)   log.Size  
#  -0.9241       1.9229  
# Random effects:
#  (Intercept)
#  1997  0.251538742
#  1998  0.803434142
#  1999 -0.091840931
#  2000  0.496705881
#  2001  0.406750132
#  2003 -0.964128413
#  2004 -0.061124277
#  2005  0.333005668
#  2006 -0.337644021
#  2007 -0.287244028
#  2008 -0.315315045
#  2009 -0.033125262
#  2010 -0.005643386
#  2011  0.110009975
#  2012 -0.247176599
log10(attributes(VarCorr(modF1))$sc)  # -0.3115011

# plots modF0
{
 # fitting
 data.f <- transform(data.f, sum.Size = log.Size+Fert)
 counts <- table(data.f$sum.Size)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 x <- seq(min(data.f$log.Size, na.rm = T), max(data.f$log.Size, na.rm = T), 0.1)
 plot(data.f$log.Size, data.f$Fert, pch = 16, col = "grey", 
  cex = sqrt(counts.col/pi), xlim = c(min(x), max(x)), main = "Fecundity", 
  ylab = "Number of flowers", xlab = "log.Size", bty = "l")
 lines(x, predict(modF0, newdata = list(log.Size = x), type = "response"), 
  lwd = 4)
  
 # residuals
 residuals <- residuals(modF0, "deviance")
 data.f <- transform(data.f, sum.res.0 = log.Size+residuals)
 counts <- table(data.f$sum.res.0)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.res.0[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modF0, type = "response")
 plot(fitted, residuals, ylim = c(-4, 4), pch = 16, bty = "n", col = "grey", 
  cex = sqrt(counts.col/pi), main = "Fecundity model residuals", 
  xlab = "scale(log10(Size))", ylab = "Deviance residuals")
 abline(0, 0)
}

# plots modF1
{
 # fitting
 data.f <- transform(data.f, sum.Size = log.Size+Fert)
 counts <- table(data.f$sum.Size)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.Size[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 plot(data.f$log.Size, data.f$Fert, bty = "n", xlim = c(min(x), max(x)), 
  main = "Fecundity", ylab = "Number of flowers", 
  xlab = "scale(log10(Size))", col = "grey", pch = 16, 
  cex = sqrt(counts.col/pi))
 ranef <- ranef(modF1)$Year
 rownames(ranef) <- as.numeric(rownames(ranef(modF1)$Year))
 legend <- data.frame(Year = unique(data.f$Year), 
  Year.col = rep(NA,length(unique(data.f$Year))), 
  Year.lty = rep(NA,length(unique(data.f$Year))))
 class(legend$Year.col) <- "character"
 for (i in 1:length(unique(data.f$Year))) {
  year.i <- unique(data.f$Year)[i]
  subset <- subset(data.f, data.f$Year == year.i)
  subset <- droplevels(subset)
  ranef.year <- ranef[which(as.numeric(rownames(ranef)) == year.i), 1]
  x <- seq(min(subset$log.Size), max(subset$log.Size), 0.1)
  y <- exp(predict(modF1, newdata = list(log.Size = x), re.form = ~ 0) 
   + ranef.year)
  legend$Year.col[i] <- as.character(subset$Year.col[1])
  legend$Year.lty[i] <- subset$Year.lty[1]
  lines(x, y, col = alpha(legend$Year.col[i], 0.5), 
   lty = legend$Year.lty[i], lwd = 2)
  #lines(x[-length(x)], y[-length(y)], col = subset$Year.col[1])
  # arrows(x[length(x)-1], y[length(y)-1], x[length(x)], y[length(y)], 
  # col = subset$Year.col[1], length = 0.1)
 }
 x <- seq(min(data.f$log.Size, na.rm = T), max(data.f$log.Size, na.rm = T), 
  0.1)
 lines(x, predict(modF1, newdata = list(log.Size = x), type = "response", 
  re.form = ~ 0), lwd = 4)  
 legend <- legend[order(legend$Year),]
 legend("topleft", legend = c("All",legend[, 1]), col = c("black", 
  alpha(legend[, 2], 0.5)), lty = c(1, legend[, 3]), bty = "n", 
  lwd = c(4, rep(2, length(legend$Year))))
    
 # residuals
 residuals <- residuals(modF1, "deviance")
 data.f <- transform(data.f, sum.res.1 = log.Size+residuals)
 counts <- table(data.f$sum.res.1)
 counts.col <- rep(1, dim(data.f)[1])
 sub.counts <- subset(counts, counts > 1)
 sum.sub.counts <- round(as.numeric(names(sub.counts)), 6)
 for (i in 1:length(sum.sub.counts))
  for (j in 1:dim(data.f)[1])
   if (round(data.f$sum.res.1[j], 6) == sum.sub.counts[i])
    counts.col[j] <- sub.counts[i]
 fitted <- predict(modF1, type = "response")
 plot(fitted, residuals, col = "grey", ylim = c(-4, 4), 
  pch = 16, cex =  sqrt(counts.col/pi), bty = "n", 
  main = "Fecundity model residuals", 
  xlab = "scale(log10(Size))", ylab = "Deviance residuals")
 abline(0, 0)
}

## Seedling size distribution analysis

# Normal on log.Size

modNS0 <- lm(log.Size ~ 1, data = data.ns)
# AIC = 1077.664
# (Intercept)
# -1.14724
log10(summary(modNS0)$sigma)  # -0.2266537

modNS1 <- lmer(log.Size ~ 1 + (1|Year), data = data.ns, REML = F)
# AIC = 979.5127
# Fixed effects:
#  (Intercept)  
#  -1.102
log10(sigma(modNS1))  # -0.2757472
# Random effects:
#  (Intercept)
#  1998  0.1734256023
#  1999 -0.0331609379
#  2000 -0.1179537375
#  2001  0.0656718643
#  2003  0.8186383185
#  2004 -0.0532607905
#  2005 -0.1663361702
#  2006 -0.2475309654
#  2007 -0.3869396839
#  2008 -0.0073740299
#  2009  0.2022492123
#  2010  0.0004822956
#  2011 -0.2523187634
#  2012  0.0044077856

# plot
x <- seq(min(data.ns$log.Size), max(data.ns$log.Size), length.out = 50)
inc <- (max(data.ns$log.Size)-min(data.ns$log.Size))/100
breaks <- c(x-inc,x[50]+inc)
hist(data.ns$log.Size, breaks = breaks, freq = F, border = F, col = "grey", 
 main = "Newborns distribution - Gaussian fitting", xlab = "scale(log10(size))", 
 ylab = "Frequency")
lines(x, dnorm(x, fixef(modNS1), sigma(modNS1)))


# Gamma on log.Size
data.ns <- transform(data.ns, log.Size.loc = log.Size - min(data.ns$log.Size) 
 + 1e-2)

modNS0 <- glm(log.Size.loc ~ 1, data = data.ns, family = Gamma)
# AIC = -296.4907
# (Intercept)  
# 2.121
summary(modNS0)$dispersion  # 1.583663

modNS1 <- glmer(log.Size.loc ~ 1 + (1|Year), data = data.ns, family = Gamma)
# AIC = -318.4779
# Fixed effects:
#  (Intercept)  
#  2.565  
# Random effects:
#  (Intercept)
#  1998 -1.11642521
#  1999 -0.48731650
#  2000 -0.01827987
#  2001 -0.82630171
#  2003 -1.84889260
#  2004 -0.34322097
#  2005  0.40282374
#  2006  1.20701583
#  2007  4.41535530
#  2008 -0.56498869
#  2009 -1.17783284
#  2010 -0.57547279
#  2011  1.29175873
#  2012 -0.50704899
attributes(VarCorr(modNS1))$sc  # 1.381432
# according to http://people.stat.sfu.ca/~raltman/stat402/402L25.pdf dispersion parameter = 1/alpha
# therefore alpha = 1/1.381432 = 0.7238864

# plot
x <- seq(min(data.ns$log.Size), max(data.ns$log.Size), length.out = 50)
inc <- (max(data.ns$log.Size)-min(data.ns$log.Size))/100
breaks <- c(x-inc,x[50]+inc)
x.loc <- x - min(data.ns$log.Size) + 1e-2
hist(data.ns$log.Size, breaks = breaks, freq = F, border = F, col = "grey", 
 main = "Newborns distribution - Gamma fitting", xlab = "scale(log10(size))", 
 ylab = "Frequency")
lines(x, dgamma(x.loc, 1/attributes(VarCorr(modNS1))$sc, 
 1/(attributes(VarCorr(modNS1))$sc)*fixef(modNS1)))

# Poisson on Size

modNS0 <- glm(Size-1 ~ 1, family = poisson, data = data.ns)
# AIC = 2119.584
# (Intercept)  
# 0.07554

modNS1 <- glmer(Size-1 ~ 1 + (1|Year), family = poisson, data = data.ns)
# AIC = 1775.972
# Fixed effects:
#  (Intercept)  
#  -0.0685 
# Random effects:
#  (Intercept)
#  1998  0.40826217
#  1999  0.05373879
#  2000 -0.23552755
#  2001  0.32879301
#  2003  1.74090757
#  2004 -0.09016732
#  2005 -0.28347872
#  2006 -0.81141198
#  2007 -1.28470906
#  2008  0.12025990
#  2009  0.53349519
#  2010 -0.03308869
#  2011 -0.38070891
#  2012  0.16006452

# plot
x <- seq(min(data.ns$Size), max(data.ns$Size))
inc <- 0.5
breaks <- c(x-inc,max(x)+inc)
x.loc <- x - 1
hist(data.ns$Size, breaks = breaks, freq = F, border = F, col = "grey",
 main = "Newborns distribution - Poisson fitting", xlab = "size", 
 ylab = "Frequency")
lines(x, dpois(x.loc, exp(fixef(modNS1))))

# Neg-binomial on Size

modNS0 <- glm.nb(Size-1 ~ 1, data = data.ns)
# AIC = 1697.736
# (Intercept)  
# 0.07554  

modNS1 <- glmer.nb(Size-1 ~ 1 + (1|Year), data = data.ns)
# AIC = 2147.872
# Fixed effects:
#  (Intercept)  
#  -0.01945  
# Random effects:
#  (Intercept)
#  1998  0.08987728
#  1999  0.00287970
#  2000 -0.08977152
#  2001  0.05863264
#  2003  0.52367698
#  2004 -0.02339985
#  2005 -0.05188238
#  2006 -0.18447376
#  2007 -0.38156393
#  2008  0.01629666
#  2009  0.12984994
#  2010 -0.01111302
#  2011 -0.07736461
#  2012  0.00840406

# plot
x <- seq(min(data.ns$Size), max(data.ns$Size))
inc <- 0.5
breaks <- c(x-inc,max(x)+inc)
x.loc <- x - 1
hist(data.ns$Size, breaks = breaks, freq = F, border = F, col = "grey",
 main = "Newborns distribution - Negative binomial fitting", 
 xlab = "size", ylab = "Frequency")
lines(x, dnbinom(x.loc, size = attributes(VarCorr(modNS1))$sc, 
 mu = exp(fixef(modNS1))))

## Flower to seedling probability
p <- dim(data.ss)[1]/sum(data.f$Fert)
# 0.09635941
log(p) + fixef(modF)[1]
# -2.903047

## Min reproductive log.Size
min.f <- min(subset(data.f$log.Size, data.f$Fert > 0))
# -1.608776

x <- seq(min(data$log.Size, na.rm = T), max(data$log.Size, na.rm = T), 
 length.out = 50)
f <- subset(data.f$log.Size, data.f$Fert > 0)
hist(f, breaks = x, freq = T, col = alpha("red", 0.5), 
 border = F, main = "Reproductive (red) vs Seedlings (blue)", 
 xlab = "log.Size")
hist(data.ns$log.Size, breaks = x, freq = T, col = alpha("blue", 0.5), 
 border = F, add = T)

x <- seq(min(data$Size, na.rm = T), max(data$Size, na.rm = T))
f <- subset(data.f$Size, data.f$Fert > 0)
hist(f, breaks = x, freq = T, col = alpha("red", 0.5),
 border = F, main = "Reproductive (red) vs Seedlings (blue)", 
 xlab = "Size")
hist(data.ns$Size, breaks = x, freq = T, col = alpha("blue", 0.5),
 border = F, add = T)