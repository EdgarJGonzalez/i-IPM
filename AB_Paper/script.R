# ANALYSIS OF THE AGE (Edad) VS BASAL AREA (AB) RELATIONSHIP

library(ggplot2)
library(gridExtra)
library(lme4)
library(mgcv)
library(reshape2)

main.dir = "/Users/Edgar/Documents/UNAM/Postdoc/Ben/i-IPM/AB_Paper"
setwd(main.dir)
data <- read.csv("datos.csv")
data$Plot <- as.factor(data$Plot)

# Plot colours
palette <- data.frame(levels(data$Plot),rainbow(length(levels(data$Plot))))
names(palette) <- c("Plot","Plot.col")
data <- merge(data, palette)

m.edad <- min(data$Edad)
M.edad <- max(data$Edad)

# Cleaning the data
data.clean = c()
for (sp in levels(data$Sp)) {
 data.sp <- subset(data, data$Sp == sp)
 for (plot in levels(data.sp$Plot))
  if (all(subset(data.sp, data.sp$Plot == plot)$AB == 0)) {
   warning(paste("plot", plot, "in sp", sp, "is empty"))
   data.sp = data.sp[-which(data.sp$Plot == plot),]
  }
 data.sp <- droplevels(data.sp)
 data.clean <- rbind(data.clean,data.sp)
}
data <- data.clean
data.pos <- subset(data, data$AB > 0)
data.pos <- droplevels(data.pos)

# MODELS

# Poly2: log10(AB>0) ~ poly(Edad,2) + (poly(Edad,2)|Plot)
for (sp in levels(data$Sp)) {
 data.sp <- subset(data, data$Sp == sp)
 data.sp <- droplevels(data.sp)
 
 data.pos.sp <- subset(data.pos, data.pos$Sp == sp)
 data.pos.sp <- droplevels(data.pos.sp)
 data.pos.sp <- transform(data.pos.sp, log10.AB = log10(AB))
 
 for (plot in levels(data.pos.sp$Plot))
  if (dim(subset(data.pos.sp, data.pos.sp$Plot == plot))[1] < 4)
   data.pos.sp <- data.pos.sp[-which(data.pos.sp$Plot == plot), ]
 if(dim(data.pos.sp)[1] == 0)
  next
 data.pos.sp <- droplevels(data.pos.sp)
 
 # models
 if (length(levels(data.pos.sp$Plot)) == 1)
  mod <- lm(log10.AB ~ poly(Edad, 2), data = data.pos.sp)
 
 if (length(levels(data.pos.sp$Plot)) > 1) {
  mod <- lmer(log10.AB ~ poly(Edad, 2) + (poly(Edad, 2) | Plot), 
   data = data.pos.sp, REML = F)
  # weigths = varIdent(form =  ~ 1 | data.pos.sp$Plot)
  if (!is.null(mod@optinfo$conv$lme4$messages))
   next	  
 }
 
 # plots
 if (length(levels(data.pos.sp$Plot)) == 1) {
  Edad.seq <- seq(min(data.pos.sp$Edad), max(data.pos.sp$Edad), 0.5)
  log10.AB.est <- coefficients(mod)[1] +
   coefficients(mod)[2]*poly(Edad.seq, 2)[, 1] +
   coefficients(mod)[3]*poly(Edad.seq, 2)[, 2]	
  data.pos.sp <- transform(data.pos.sp, poly2.res = residuals(mod), 
   type = "obs") 
  new.data <- as.data.frame(cbind(Edad.seq, log10.AB.est, "est", 0, 
   "black"))
  data.plot <- data.frame(data.pos.sp$Edad, data.pos.sp$log10.AB, 
   data.pos.sp$type, data.pos.sp$Plot, data.pos.sp$Plot.col)
  names(data.plot) <- names(new.data) <- c("Edad", "log10.AB", "type", "Plot", 
   "Plot.col")
  data.plot <- rbind(data.plot, new.data) 
  data.plot <- transform(data.plot,Edad = as.numeric(Edad),
   log10.AB = as.numeric(log10.AB))
  
  plot1 <- ggplot(data.plot,aes(Edad,log10.AB)) + 
  geom_line(data = data.plot[data.plot$type == "obs", ], aes(group = Plot),
   colour = data.plot[data.plot$type == "obs", ]$Plot.col) +
  geom_line(data = data.plot[data.plot$type == "est", ], aes(group = Plot), 
   colour = "black") +
  theme(panel.background = element_rect(fill = "white")) + 
  ggtitle(sp) +
  scale_x_continuous(limits = c(m.edad, M.edad)) +
  ylab("log10 Basal area (units)") +
  xlab("Age (year)")
   
  plot2 <- ggplot(data.pos.sp, aes(Plot,poly2.res)) + 
  geom_boxplot(colour = levels(data.pos.sp$Plot.col)) +
  theme(panel.background = element_rect(fill = "white")) + 
  ylab("Poly2 residuals") +
  xlab("Plot")
 }
 
 if (length(levels(data.pos.sp$Plot)) > 1) {
  Edad.seq <- seq(min(data.pos.sp$Edad),max(data.pos.sp$Edad), 0.5)
  log10.AB.est <- fixef(mod)[1] + fixef(mod)[2]*poly(Edad.seq, 2)[, 1] +
   fixef(mod)[3]*poly(Edad.seq, 2)[, 2]
  data.pos.sp <- transform(data.pos.sp, poly2.res = residuals(mod), 
   type = "obs")
  new.data <- as.data.frame(cbind(Edad.seq, log10.AB.est, "fix.est", 0, 
   "black"))
  data.plot <- data.frame(data.pos.sp$Edad, data.pos.sp$log10.AB, 
   data.pos.sp$type, data.pos.sp$Plot, data.pos.sp$Plot.col)
  names(data.plot) <- names(new.data) <- c("Edad", "log10.AB", "type", "Plot", 
   "Plot.col")
  data.plot <- rbind(data.plot, new.data)
  new.data = c()
  for(plot in levels(data.pos.sp$Plot)) {
   data.pos.sp.plot <- subset(data.pos.sp, data.pos.sp$Plot == plot)
   data.pos.sp.plot <- droplevels(data.pos.sp.plot)
   Edad.seq <- seq(min(data.pos.sp.plot$Edad), max(data.pos.sp.plot$Edad), 
    0.25)
   coefs <- as.numeric(subset(coef(mod)$Plot, 
    rownames(coef(mod)$Plot) == plot))
   log10.AB.plot.est <- coefs[1] + coefs[2]*poly(Edad.seq, 2)[, 1] + 
    coefs[3]*poly(Edad.seq, 2)[, 2]  # coefs includes fix and ran effects
   new.data <- rbind(new.data, data.frame(Edad.seq, log10.AB.plot.est, 
	"ran.est", plot, levels(data.pos.sp.plot$Plot.col)))
  }
  names(new.data) <- c("Edad", "log10.AB", "type", "Plot", 
   "Plot.col")
  data.plot <- rbind(data.plot, new.data)
  data.plot <- transform(data.plot,Edad = as.numeric(Edad),
   log10.AB = as.numeric(log10.AB))
  
  plot1 <- ggplot(data.plot, aes(Edad,log10.AB)) + 
  geom_line(data = data.plot[data.plot$type == "obs", ], aes(group = Plot), 
   colour = data.plot[data.plot$type == "obs", ]$Plot.col) +
  geom_line(data = data.plot[data.plot$type == "ran.est", ], aes(group = Plot), 
   colour = data.plot[data.plot$type == "ran.est", ]$Plot.col, alpha = 0.25) +
  geom_line(data = data.plot[data.plot$type == "fix.est", ], colour = "black") +
  theme(panel.background = element_rect(fill = "white")) + 
  ggtitle(sp) +
  scale_x_continuous(limits = c(m.edad, M.edad)) +
  ylab("log10 Basal area (units)") +
  xlab("Age (years)")
  
  plot2 <- ggplot(data.pos.sp, aes(Plot,poly2.res)) + 
  geom_boxplot(aes(group = Plot), colour = levels(data.pos.sp$Plot.col)) +
  theme(panel.background = element_rect(fill = "white")) + 
  ylab("Poly2 residuals") +
  xlab("Plot")
 }
 
 pdf(paste("poly2_poly2_", sp, ".pdf", sep = ""))
 grid.arrange(plot1, plot2, nrow = 2)
 dev.off()
 
}
