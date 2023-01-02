#######################
###### LIBRARIES ######
#######################

#this is a test comment for git
#test comment number 2

library(rstudioapi)
library(dplyr)
library(lubridate)
library(lme4)
library(ggplot2)
library(ggeffects)
library(DHARMa)
library(MuMIn)
library(nlstools)

# Look in script folder for data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(na.action="na.fail")


#######################
###### 1. GROWTH ######
#######################

### 1a. Fitting Gompertz growth curve to age vs size data ###

# Define the Gompertz function
gz2 <- size.mm~larvalsize.mm*exp(log(adult.size.mm/larvalsize.mm)*(1-exp(-kg*age)))

# Read data
fa <- read.csv("size_age.csv")

# Subset the data to females only and filter out NAs, separate treatments for convenience
fha <- subset(fa,sex=='f' & treatment=='h' & !is.na(size.mm) & !is.na(age))
fla <- subset(fa,sex=='f' & treatment=='l' & !is.na(size.mm) & !is.na(age))

# Non-linear least squares fit for Gompertz parameters
# For high-food
gstarts <- c(larvalsize.mm = 0.08, kg = 0.23, adult.size.mm = 0.7)
gfit.h <- nls(gz2,data=fha,start=gstarts)
overview(gfit.h)

# For low-food
gfit.l <- nls(gz2,data=fla,start=gstarts)
overview(gfit.l)

# Compare the models
anova(gfit.h,gfit.l)

# Extract and model the variance
# H
gvh <- nlsResiduals(gfit.h)[["resi1"]]
gvh[,2] <- gvh[,2]^2
plot(gvh[,1],gvh[,2])
hvar <- gvh[,2]; hsizes<-gvh[,1]
v01 <- lm(hvar~hsizes+0); summary(v01)

# L
gvl <- nlsResiduals(gfit.l)[["resi1"]]
gvl[,2] <- gvl[,2]^2
plot(gvl[,1],gvl[,2])
lvar <- gvl[,2]; lsizes<-gvl[,1]
v11 <- lm(lvar~lsizes+0); summary(v11)

# Plotting
ages <- 0:70

# Parameters extracted
# H
w0.h <- 0.08631
kg.h <- 0.25111
A.h <- 0.69896
gstarts.h <- c(w0.h = 0.1, kg.h = 0.2, A.h = 0.7)

#L
w0.l <- 0.08815
kg.l <- 0.21725
A.l <- 0.72186
gstarts.l <- c(w0.l = 0.1, kg.l = 0.2, A.l = 0.7)

# Plot together and save
# H is red, L is blue
ggplot(fla,aes(x=age,y=size.mm))+
  geom_point(color="#377EB8",alpha=0.6)+
  geom_smooth(method="nls",
              formula=y~A.l*exp((log(w0.l/A.l))*(exp(-kg.l*x))),
              method.args=list(start=gstarts.l),
              se=FALSE,color="#377EB8") +
  geom_point(data=fha,aes(x=age,y=size.mm),color="#E41A1C",alpha=0.6)+
  geom_smooth(data=fha,aes(x=age,y=size.mm),method="nls",
              formula=y~A.h*exp((log(w0.h/A.h))*(exp(-kg.h*x))),
              method.args=list(start=gstarts.h),
              se=FALSE,color="#E41A1C") +
  labs(x="\nAge (days)",y="Size (mm)\n",tag="a.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",legend.box.background = element_rect(color="black", size=1))
ggsave("gompertz.png", dpi=300, height=3.1, width=3.1, units="in")


### 1b. Plotting size_t vs size_t+1 ###

# Read in exact fits and log them
g <- read.csv("gomp_lin.csv")
g$lt <- log(g$t); g$lt1 <- log(g$t1)

# Plot and save
ggplot(data=g,aes(x=lt,y=lt1,color=Treatment)) +
  geom_smooth(method="lm",se=FALSE,fullrange=FALSE) +
  labs(x="\nlog Size (t)",y="log Size (t+1)\n") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",legend.box.background = element_rect(color="black", size=1)) +
  scale_colour_brewer(palette="Set1")
ggsave("log_size.png", dpi=300, height=4, width=4, units="in")


#########################
###### 2. SURVIVAL ######
#########################

# Read in daily time-step data with interpolated growth, set categorical variables to factors, remove rows with missing survival data
mats <- read.csv("ipm_gompertz.csv")
mats$survival <- as.factor(mats$survival); mats$treatment <- as.factor(mats$treatment); mats$concat <- as.factor(mats$concat)
mats_s <- subset(mats,!is.na(survival))

# Specify the full model and perform full-model evaluation based on AIC using dredge()
# full.s <- glmer(survival ~ size*treatment+census*treatment +I(census^2)*treatment + (1|concat),family=binomial,nAGQ=0,data=mats_s);summary(full.s)
# dredge(full.s)

# Compare top models
final.s <- glmer(survival ~ size+census*treatment+I(census^2) + (1|concat),family=binomial,nAGQ=0,data=mats_s);summary(final.s)
final.s2 <- glmer(survival ~ size+census+I(census^2) + (1|concat),family=binomial,nAGQ=0,data=mats_s);summary(final.s2)
anova(final.s,final.s2)

# Plot and save
mats_s$sv <- as.double(mats_s$survival) -1
ggpredict(final.s2, c("size[all]", "census[5,15,30,60]"), ci.lvl=NA) |> plot(add.data = FALSE, limit.range = TRUE) +
  scale_y_continuous(breaks = seq(.7, 1.05, by = .1)) +
  geom_boxplot(data=mats_s,width=0.02,aes(x=size,y=ifelse(sv>0,sv+.05,sv+.7),group=sv),inherit.aes=FALSE) +
  labs(x="\nSize (mm)",y="Survival\n",color="Age",title=NULL,tag="b.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",legend.box.background = element_rect(color="black", size=1),
        axis.line.x.bottom=element_line(color="black"),axis.line.y.left=element_line(color="black"),
        axis.ticks=element_line(color="black"),axis.text = element_text(face="bold"),
        axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black"))
ggsave("survival.png", dpi=300, height=3.1, width=3.1, units="in")


#############################
###### 3. REPRODUCTION ######
#############################

# Incompatible with preceding models
library(lmerTest)

### 3a. Plotting inter-clutch time for supps ###

# Get individual level data
f <- read.csv("vr_individuals.csv")
r <- read.csv("vr_reproduction.csv")

# Map chemostat/population number to treatment
highs <- c(1, 4, 5, 6, 11, 12, 13, 15, 17, 18)
lows <- c(2, 7, 8, 9, 10, 14, 16, 19, 20)

# Create primary key for matching across tables
f$concat <- paste(f$lineage,f$focal.ID,f$round,sep=".")
r$concat <- paste(r$Lineage,r$Mother.ID,r$Round,sep=".")

# Calculate mother age for each clutch
r$m.t0 <- f$collection.date[match(r$concat, f$concat)]
r$t.hatch <- as.numeric(dmy(r$Hatch.Date)-dmy(r$m.t0))

# Fill in treatment and flag females as whether they ever reproduced (or not)
f <- f %>%
  mutate(treatment = case_when(
    lineage %in% highs ~ 'h',
    lineage %in% lows ~ 'l'
  )) %>%
  mutate(reproducer = case_when(
    f$concat %in% r$concat ~ 1,
    TRUE ~ 0
  ))

r <- r %>%
  mutate(treatment = case_when(
    Lineage %in% highs ~ 'h',
    Lineage %in% lows ~ 'l'
  ))

# Infill final mother size
r$fin.sz <- f$adult.size.mm[match(r$concat, f$concat)]

# Filter out rows with missing data for: hatch date, offspring number, brood time, mother age and mother size
r.brood <- subset(r,!is.na(Hatch.Date) & !is.na(Offspring.Count))
r.brood <- subset(r.brood,Offspring.Count > 0)
r.brood$brood.time <- as.numeric(dmy(r.brood$Hatch.Date)-dmy(r.brood$Collection.Date))
r.brood <- subset(r.brood,!is.na(brood.time) & brood.time >= 0)
r.brood <- subset(r.brood,!is.na(fin.sz) & !is.na(t.hatch))

# Plot and save
ggplot(r.brood,aes(x=brood.time,fill=treatment)) +
  geom_histogram(alpha=0.6,binwidth=1,position='identity') +
  labs(x="\nBrood time (days)",y="Count\n") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 80, by = 10)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",legend.box.background = element_rect(color="black", size=1)) +
  scale_fill_brewer(palette="Set1")
ggsave("brood_time.png", dpi=300, height=4, width=4, units="in") #72dpi and 300x300=4.16 x 4.16 in


### 3b. Probability of reproducing ###

# Subset daily timestep data (from SURVIVAL) to reproductive females only, encode aborted clutches as binary (for later)
mat <- mats[which(mats$sex=='f' | mats$concat %in% r$concat),]
mat <- mat %>%
  mutate(blanks = case_when(
    reproduction == 1 & fecundity == 0 ~ TRUE,
    reproduction == 1 & fecundity > 0 ~ FALSE,
    TRUE ~ NA)
  )

# Filter rows after death
mat_b <- subset(mat,!is.na(survival))

# Full model and selecion based on AIC using dredge()
#full.r <- glmer(reproduction ~ treatment*size+treatment*census+I(census^2)*treatment + (1|concat),family=binomial,nAGQ=0,na.action="na.fail",data=mat_b); summary(full.r)
#dredge(full.r)

# Best model
final.nu <- glmer(reproduction ~ size+treatment*census+treatment*I(census^2) + (1|concat),family=binomial,nAGQ=0,na.action="na.omit",data=mat_b); summary(final.nu)

# Plot and save
pred.r <- subset(mat_b,!is.na(reproduction))
pred.r$pred.r <- predict(final.nu,re.form=NA)
pred.r$treatment <- as.factor(pred.r$treatment)
levels(pred.r$treatment) <- c("High food","Low food")

# L panel
ggpredict(final.nu, c("size[all]", "census[18,23,28,33]", "treatment[l]"), ci.lvl=NA) |> plot(add.data =FALSE, limit.range = TRUE) +
  scale_y_continuous(breaks = seq(0, .25, by = .05)) +
  geom_boxplot(data=subset(pred.r,treatment=="Low food"),width=0.02,aes(x=size,y=ifelse(reproduction>0,reproduction-.7,reproduction-.05),group=reproduction),inherit.aes=FALSE) +
  labs(x="\nMaternal size (mm)",y="Probability of producing eggs\n",color="Age",title=NULL,tag="c.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none", legend.box.background = element_rect(color="black", size=1),
        axis.line.x.bottom=element_line(color="black"),axis.line.y.left=element_line(color="black"),
        axis.ticks=element_line(color="black"),axis.text = element_text(face="bold"),
        axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black"))
ggsave("prob clutch l2.png", dpi=300, height=3.1, width=3.1, units="in")

# H panel
ggpredict(final.nu, c("size[all]", "census[18,23,28,33]", "treatment[h]"), ci.lvl=NA) |> plot(add.data =FALSE, limit.range = TRUE) +
  scale_y_continuous(breaks = seq(0, .25, by = .05)) +
  geom_boxplot(data=subset(pred.r,treatment=="High food"),width=0.02,aes(x=size,y=ifelse(reproduction>0,reproduction-.7,reproduction-.05),group=reproduction),inherit.aes=FALSE) +
  labs(x="\nMaternal size (mm)",y="Probability of producing eggs\n",color="Age",title=NULL,tag="d.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none", legend.box.background = element_rect(color="black", size=1),
        axis.line.x.bottom=element_line(color="black"),axis.line.y.left=element_line(color="black"),
        axis.ticks=element_line(color="black"),axis.text = element_text(face="bold"),
        axis.title.x=element_text(colour="black"), axis.title.y=element_text(colour="black"))
ggsave("prob clutch h2.png", dpi=300, height=3.1, width=3.1, units="in")


### 3c. Probability of (not) aborting a clutch ###

# Filter rows to living mums without missing data on clutch abortion
mat$blanks <- as.factor(mat$blanks)
mat_c <- subset(mat,!is.na(survival) & !is.na(blanks))

# Full model and model selection based on AIC with dredge()
#br <- glmer(blanks ~ size*census*treatment+I(census^2)*treatment + (1|concat),family=binomial,nAGQ=0,na.action="na.fail",data=mat_c); summary(br)
#dredge(br)

# Final model
final.b2 <- glmer(blanks ~ size + (1|concat),family=binomial,nAGQ=0,na.action="na.fail",data=mat_c); summary(final.b2)

# Plot and save
pred.b <- mat_c #subset(mat,!is.na(blanks))
pred.b$blanks <- as.integer(pred.b$blanks) -1
pred.b$pred.b <- predict(final.b2,re.form=NA)

ggplot(data=pred.b,aes(x=size,y=blanks)) +
  stat_smooth(method="glm",se=FALSE,fullrange=FALSE,
              method.args = list(family=binomial)) +
  geom_boxplot(data=pred.b,width=0.03,aes(x=size,y=ifelse(blanks>0,blanks-.5,blanks-.05),group=blanks)) +
  labs(x="\nMaternal size (mm)",y="Probability of egg failure\n",color="Age",tag="e.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.box.background = element_rect(color="black", size=1))
ggsave("recruitment.png", dpi=300, height=3.1, width=3.1, units="in")


### 3d. Clutch size ###

# From individual reproduction data, filter out rows with missing entries for cohort number, offspring count, hatch date, mother size
r2 <- r[!is.na(r$Cohort.num),]; r2 <- r2[!is.na(r2$Offspring.Count),]; r2 <- r2[!is.na(r2$t.hatch),]; r2 <- r2[!is.na(r2$fin.sz),] 

# Only retain clutches that yielded offspring
r2z <- subset(r2,Offspring.Count>0)

# Full model with AIC full model selection with dredge()
#full.n <- glmer(Offspring.Count ~ fin.sz*treatment + t.hatch*treatment + I(t.hatch^2)*treatment + (1|concat),family=poisson,nAGQ=0,na.action="na.fail",data=r2z); summary(full.n)
#dredge(full.n)

# Final model
final.n <- glmer(Offspring.Count ~ t.hatch*treatment+I(t.hatch^2) + (1|concat),family=poisson,nAGQ=0,data=r2z); summary(final.n)

# Plot and save
pred.f <- subset(r2z,!is.na(Offspring.Count))
ldat <- data.frame(t.hatch=15:45,treatment="l"); hdat <- data.frame(t.hatch=15:45,treatment="h"); ndat <- rbind(ldat,hdat)
ndat$pred <- predict(final.n,newdata=ndat,type="response",re.form=NA)

ggplot(data=pred.f,aes(x=t.hatch,y=Offspring.Count,color=treatment)) +
  geom_line(data=ndat,aes(y=pred,x=t.hatch,color=treatment)) +
  geom_point(alpha = 0.6) + 
  labs(x="\nMaternal age (days)",y="Clutch size\n",tag="f.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",legend.box.background = element_rect(color="black", size=1)) +
  scale_colour_brewer(palette="Set1")
ggsave("fecundity.png", dpi=300, height=3.1, width=3.1, units="in")


### 3e. Egg size ###

# Aggregate reproductive data to clutch
r.ag <- r2 %>%
  group_by(concat,Cohort.num,treatment,Offspring.Count,fin.sz,t.hatch,Eggsac.area.mm2) %>%
  summarise(eg.av=mean(c(egg1,egg2,egg3,egg4,egg5,egg6,egg7,egg8,egg9,egg10),na.rm=TRUE),n.av=mean(c(larva1,larva2,larva3,larva4,larva5,larva6,larva7,larva8,larva9,larva10),na.rm=TRUE))

#Remove clutches that produced no offspring
rag.z <- subset(r.ag,Offspring.Count > 0)

# Full model with AIC full model selection with dredge()
#full.e <- lmer(eg.av~ treatment*fin.sz*t.hatch+treatment*I(t.hatch^2)+(1|concat), REML=FALSE,na.action="na.fail",data=rag.z); summary(full.e)
#dredge(full.e)

# LRT evaluation of best models
final.e <- lmerTest::lmer(eg.av~ as.factor(treatment)+(1|concat), REML=FALSE, data=rag.z); summary(final.e)
final.e2 <- lmerTest::lmer(eg.av~ 1+(1|concat), REML=FALSE, data=rag.z); summary(final.e2)

anova(final.e,final.e2)

#Modelling variance
gv <- resid(final.e)^2
treatment <- final.e@frame$'as.factor(treatment)'
plot(treatment,gv)

v0 <- lm(gv~1); summary(v0)

# Bootstrap 95% confidence intervals
bb <- bootMer(final.e,
              FUN=function(x)predict(x, re.form=NA),
              nsim=500)

pred.e <- rag.z
pred.e$lci <- apply(bb$t, 2, quantile, 0.025)
pred.e$uci <- apply(bb$t, 2, quantile, 0.975)
pred.e$pred.e <- predict(final.e,re.form=NA) 

pred.e %>%
  group_by(treatment) %>%
  summarise(mean = mean(pred.e),
            hci = mean(uci),
            dci = mean(lci))

# Plot and save
pred.e$treatment <- as.factor(pred.e$treatment)
levels(pred.e$treatment) <- c("High food","Low food")
ggplot(pred.e, aes(x=treatment, y=pred.e, group=treatment)) + 
  geom_point(stat="identity", position=position_dodge(0.1), aes(colour=treatment)) +
  geom_errorbar(aes(ymin=lci, ymax=uci, colour=treatment), width=.2, position=position_dodge(.1)) +
  labs(x="\nLineage",y="Mean egg diameter (mm)\n",tag="g.") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",legend.box.background = element_rect(color="black", size=1)) +
  scale_colour_brewer(palette="Set1")
ggsave("inheritance.png", dpi=300, height=3.1, width=3.1, units="in")