# FITTING AN IPM TO Roberto Salguero-Gomez DATA

require(lme4)
require(scales)

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
    z.Size = ifelse(is.na(Size), NA, scale(log10(Size), scale = sd.Size)))


## Survival analysis
data.s <- c()
for (ind in levels(data$ID)) {
    data.ind <- subset(data, data$ID == ind & is.na(data$z.Size) != T)
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

modS0 <- glmer(Survival ~ poly(z.Size,2) + (1|Plot/Quadrat) + (1|Year), 
    family = binomial, data = data.s)
# AIC = 3448.568

modS1 <- glmer(Survival ~ poly(z.Size,2) + (1+z.Size|Plot/Quadrat) + 
    (1+z.Size|Year), family = binomial, data = data.s)
# Model failed to converge with max|grad| = 0.00254122 (tol = 0.001, component 11)
# AIC = 3433.626
# Fixed Effects:
#  (Intercept)  poly(z.Size, 2)1  poly(z.Size, 2)2  
#  0.2301       28.3015           -2.5211  

modS2 <- glmer(Survival ~ poly(z.Size,2) + (poly(z.Size,2)|Plot/Quadrat) + 
    (poly(z.Size,2)|Year), family = binomial, data = data.s)
# AIC = 3435.640

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
plot(data.s$z.Size, data.s$Survival, pch = 16, col = alpha("blue", 0.25), 
    main = "Survival", ylab = "Survival probability", xlab = "z.Size", 
    bty = "l")
lines(x, predict(modS1, newdata = list(z.Size = x), type = "response", 
    re.form = ~ 0))

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
y <- predict(modS1, newdata = list(z.Size = x), type = "response", 
    re.form = ~ 0)
x <- 10^(x*sd.Size+mu.Size)
plot(data.s$Size, data.s$Survival, pch = 16, col = alpha("blue", 0.25), 
    main = "Survival", ylab = "Survival probability", xlab = "Size", bty = "l")
lines(x, y)


## Growth analysis
data.G <- subset(data, is.na(data$z.Size) == F)
data.g <- c()
for (ind in levels(data.G$ID)) {
    data.ind <- subset(data.G, data.G$ID == ind)
    if (dim(data.ind)[1] > 1) {
        size1 <- data.ind$z.Size[1:(length(data.ind$z.Size)-1)]
        size2 <- data.ind$z.Size[2:length(data.ind$z.Size)]
        data.ind <- data.ind[-1,]
        data.ind <- transform(data.ind, z.Size1 = size1, z.Size2 = size2)
        data.g <- rbind(data.g, data.ind)
    }
}
data.g <- subset(data.g, data.g$Year != 2003)

modG <- lmer(z.Size2 ~ z.Size1 + (1|Plot/Quadrat) + (1|Year), data = data.g, 
    REML = F)
# AIC = 4273.905

modG <- lmer(z.Size2 ~ z.Size1 + (1|Plot/Quadrat) + (1+z.Size1|Year), 
    data = data.g, REML = F) 
# AIC = 4254.191

modG <- lmer(z.Size2 ~ z.Size1 + (1|Plot/Quadrat), data = data.g, REML = F)
# AIC = 4488.523

modG <- lmer(z.Size2 ~ z.Size1 + (1|Year), data = data.g, REML = F)
# AIC = 4309.704

modG <- lmer(z.Size2 ~ z.Size1 + (1+z.Size1|Plot/Quadrat) + (1+z.Size1|Year), 
    data = data.g, REML = F) 
# AIC = 4220.726
# Fixed Effects:
#  (Intercept)  z.Size1  
#  0.3160       0.6332
log10(attributes(VarCorr(modG))$sc)
# -0.2315721

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
plot(data.g$z.Size1, data.g$z.Size2, pch = 16, col = alpha("blue", 0.25),
    main = "Growth", ylab = "z.Size at time t+1", xlab = "z.Size at time t", 
    bty = "l" )
lines(x, predict(modG, newdata = list(z.Size1 = x), type = "response", 
    re.form = ~ 0))

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
y <- predict(modG, newdata = list(z.Size1 = x), 
    type = "response", re.form = ~ 0)
x <- 10^(x*sd.Size+mu.Size)
y <- 10^(y*sd.Size+mu.Size)
x.obs <- 10^(data.g$z.Size1*sd.Size+mu.Size)
y.obs <- 10^(data.g$z.Size2*sd.Size+mu.Size)
plot(x.obs, y.obs, pch = 16, col = alpha("blue", 0.25), main = "Growth", 
    ylab = "Size at time t+1", xlab = "Size at time t", bty = "l")
lines(x, y)


## Fecundity
data.f <- subset(data, is.na(data$Fert) == F)

modF <- glmer(Fert ~ z.Size + (1|Plot/Quadrat) + (1|Year),  family = poisson, 
    data = data.f)
# AIC = 10226.468
# Fixed Effects:
#  (Intercept)   z.Size  
#  -0.5634       1.4079

modF <- glmer(Fert ~ z.Size + (1+z.Size|Plot/Quadrat) + (1+z.Size|Year),
    family = poisson, data = data.f)
# AIC = 10240.20

modF <- glmer.nb(Fert ~ z.Size + (1+z.Size|Plot/Quadrat) + (1+z.Size|Year), data = data.f)
# AIC = 9142.683

modF <- glmer.nb(Fert ~ z.Size + (1|Plot/Quadrat) + (1|Year), data = data.f)
# AIC = 9139.437
# Fixed Effects:
#  (Intercept)   z.Size  
#  -0.7636       1.6242  
log10(fixef(modF)[2]) 
# 0.2106348

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
plot(data.f$z.Size, data.f$Fert, pch = 16, col = alpha("blue", 0.25),
    main = "Fecundity", ylab = "Number of flowers", xlab = "z.Size", bty = "l")
lines(x, predict(modF, newdata = list(z.Size = x), type = "response", 
    re.form = ~ 0))

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), 0.1)
y <- predict(modF, newdata = list(z.Size = x), 
    type = "response", re.form = ~ 0)
x <- 10^(x*sd.Size+mu.Size)
plot(data.f$Size, data.f$Fert, pch = 16, col = alpha("blue", 0.25),
    main = "Fecundity", ylab = "Number of flowers", xlab = "Size", bty = "l")
lines(x, y)


## Seedling size distribution
data.ns <- c()
for(ind in levels(data$ID)) {
    data.ind <- subset(data, data$ID == ind)
    if(dim(data.ind)[1] > 1)
        for(i in 2:length(unique(data.ind$Year)))
            if(is.na(data.ind$z.Size[i-1]) & !is.na(data.ind$z.Size[i]))
                data.ns <- rbind(data.ns, data.ind[i, ])
}
# individuals that are already reproducing on the first year they're recorded 
#  are not considered seedlings
data.ns <- subset(data.ns, data.ns$Fert == 0)

### Normal on z.Size

modNS <- lmer(z.Size ~ 1 + (1|Plot/Quadrat) + (1|Year), data = data.ns, 
    REML = F)
# AIC = 1056.0156
# Fixed Effects:
#  (Intercept)  
#  -1.094 
log10(attributes(VarCorr(modNS))$sc)
#  -0.288406

modNS <- lmer(z.Size ~ 1 + (1+z.Size|Plot/Quadrat) + (1+z.Size|Year),
    data = data.ns, REML = F)
# AIC = -9673.734
# Fixed Effects:
#  (Intercept)  
#  0.1612  
log10(attributes(VarCorr(modNS))$sc)
#  7.528599e-05
# THE MODEL COULDN'T BE CORRECTLY ESTIMATED

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T),
    length.out = 50)
hist(data.ns$z.Size, breaks = x, freq = F, border = "blue", 
    main = "Seedling size", xlab = "z.Size", ylab = "Frequency")
lines(x, dnorm(x, fixef(modNS), attributes(VarCorr(modNS))$sc))

### Poisson on Size

modNS <- glmer(Size ~ 1 + (1|Plot/Quadrat) + (1|Year), data = data.ns, 
    family = poisson)
# AIC = 1991.1533

modNS1 <- glmer(Size ~ 1 + (1|Plot/Quadrat) + (1+Size|Year), data = data.ns, 
    family = poisson)
# AIC = 1597.9609
# Fixed Effects:
#  (Intercept)
#  1.181
log10(fixef(modNS1))
# 0.07234227

x <- seq(min(data.ns$Size, na.rm = T), max(data.ns$Size, na.rm = T), 1)
y <- dpois(x, fixef(modNS1))
hist(data.ns$Size, breaks = seq(0.5,max(data.ns$Size)+0.5,1), freq = F, 
    border = "blue", main = "Seedling size", xlab = "Size", 
    ylab = "Frequency")
lines(x, y) 

## Gamma on z.Size
data.ns <- transform(data.ns, z.Size.loc = z.Size - min(data.ns$z.Size) + 1e-2)

modNS <- glmer(z.Size.loc ~ 1 + (1|Plot/Quadrat) + (1|Year), data = data.ns, 
    family = Gamma)
# Spatial structure couldn't be estimated in the model

modNS <- glmer(z.Size.loc ~ 1 + (1|Year), data = data.ns, family = Gamma)
# AIC = -3826.858
# Fixed Effects:
#  (Intercept)  
#  2.584

x <- seq(min(data.ns$z.Size, na.rm = T), max(data.ns$z.Size, na.rm = T),
        length.out = 50)
inc <- (max(data.ns$z.Size)-min(data.ns$z.Size))/100
breaks <- c(x-inc,x[50]+inc)
x.loc <- x - min(data.ns$z.Size) + 1e-2
hist(data.ns$z.Size, breaks = breaks, freq = F, border = "blue")
lines(x, dgamma(x.loc, 1/fixef(modNS)))


## Flower to seedling probability
p <- dim(data.ss)[1]/sum(data.f$Fert)
# 0.09635941
log(p) + fixef(modF)[1]
# -2.903047

## Min reproductive z.Size
min.f <- min(subset(data.f$z.Size, data.f$Fert > 0))
# -1.608776
max.ns <- max(data.ns$z.Size)

x <- seq(min(data$z.Size, na.rm = T), max(data$z.Size, na.rm = T), length.out = 50)
hist(data.f$z.Size, breaks = length(x), freq = T, border = "red",
    main = "Reproductive (red) vs Seedlings (blue)", xlab = "z.Size")
hist(data.ns$z.Size, breaks = length(x), freq = T, border = "blue", add = T)

x <- seq(min(data$Size, na.rm = T), max(data$Size, na.rm = T), 1)
hist(data.f$Size, breaks = length(x), freq = T, border = "red",
    main = "Reproductive (red) vs Seedlings (blue)", xlab = "Size")
hist(data.ns$Size, breaks = length(x), freq = T, border = "blue", add = T)
