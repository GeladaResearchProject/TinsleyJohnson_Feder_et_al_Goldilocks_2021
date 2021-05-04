## LOAD PACKAGES NEEDED
require(plyr);require(dplyr);require(lme4); require(lmerTest); require(MuMIn); require(ggplot2);require(ggplot2);require(ggthemes);require(ggbeeswarm);require(effects);require(coxme);require(reshape2)
theme_set(theme_minimal(base_size=20))

## FUNCTION FOR ASSESSING VIFs FROM LME4 OBJECTS
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# There are 8 R objects:
# dob_ranges: minimum and maximum estimated dobs; used to randomize dates within this range for survival analysis

# female_survival: monthly female data, including the number of unit females, female age, and # of males, 
# indicating whether the female survived (0) or died (1)

# female_repro_1mos: monthly female data, 
# including unit females, female age, and # of males, 
# indicating whether the female successfully reproduced (i.e., conceived a surviving infant)

# female_repro_6mos: female data grouped into six-month bins (i.e., Jan-Jun, Jul-Dec), 
# including unit females, female age, and # of males, 
# indicating whether the female successfully reproduced (i.e., conceived a surviving infant)

# ibi_data: interbirth intervals, including maternal age, # females at birth, infant sex,
# and whether a takeover occurred

# infant_survival: infant survival data, including # females at birth, infant sex, 
# maternal ID; infant end ages set at age at death, age at end of observation, or 1.5 (if known to have survived)
# "Censored" indicates whether the infant died during the interval

# unit_takeovers: unit periods of variable length
# number of females and males indicated, as well as number of months unit had a given composition
# number of takeovers indicated

# infant_mortality_causes: all infants separated by potential cause of death
# data structured similarly to infant_survival

#############################################################################
###                         ADULT FEMALE SURVIVAL                         ###
#############################################################################
female_mortality$unit.cat<-"SMALL"
female_mortality$unit.cat[female_mortality$N.Females>4.5]<-"MEDIUM"
female_mortality$unit.cat[female_mortality$N.Females>7.5]<-"LARGE"
female_mortality$unit.cat<-factor(female_mortality$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))

summary(lm(female_mortality$Age~female_mortality$N.Females))

ks.test(female_mortality$Age[female_mortality$unit.cat=="SMALL"],female_mortality$Age[female_mortality$unit.cat=="MEDIUM"])
ks.test(female_mortality$Age[female_mortality$unit.cat=="SMALL"],female_mortality$Age[female_mortality$unit.cat=="LARGE"])
ks.test(female_mortality$Age[female_mortality$unit.cat=="MEDIUM"],female_mortality$Age[female_mortality$unit.cat=="LARGE"])

library(ggridges)
ggplot(data=female_mortality, aes(x=Age, y=unit.cat, fill=unit.cat)) + stat_density_ridges(quantile_lines = T, quantile_fun = median) + theme(legend.position = "none") + scale_fill_manual(values=c("#993404","#ec7014","#fec44f")) + ylab("unit size category") + xlab("age (years)")

# Construct model
mortality.model<-glmer(Censored ~ poly(N.Females,2) + poly(Age,2) + N.Males + (1|Individual) +(1|Unit) + (1|Year), family="binomial", data=female_mortality, control=glmerControl(optimizer="bobyqa")) 
summary(mortality.model)

## MAKE MEDIAN SURVIVAL ESTIMATES
female_mortality$Age<-round(female_mortality$Age,1)
female_mortality$End.Age<-round(female_mortality$End.Age,1)
female_mortality[female_mortality$Age==female_mortality$End.Age,]
female_mortality$End.Age[female_mortality$Age==female_mortality$End.Age]<-female_mortality$End.Age[female_mortality$Age==female_mortality$End.Age]+0.01
female_mortality[female_mortality$Age>female_mortality$End.Age,]
surv_object<-Surv(time=female_mortality$Age, time2=female_mortality$End.Age, event=female_mortality$Censored)
f1 <- survfit(surv_object ~ unit.cat, data = female_mortality)
f1
f.null <- survfit(surv_object ~ 1, data = female_mortality)
f.null
plot(f.null)
plot(f1)

# Figure S2 - Survival curve across unit sizes
library(survminer)
survival2<-ggsurvplot(f1, data = female_mortality, xlim=c(5,30), ylim=c(0,1), legend="none", conf.int = F,  censor=F, xlab="age (years)", ylab="proportion adults surviving", legend.labs=c("SMALL", "MEDIUM", "LARGE"), legend.title="", break.time.by=5, palette = c("#993404","#ec7014","#fec44f"), font.x=c(size=20), font.y=c(size=20))
plot2C<-survival2$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot2C 

## RUN ANALYSIS ACROSS 1000 SIMULATIONS OF FEMALE AGE
dob.ranges[,2:4]<-lapply(dob.ranges[,2:4], as.Date, format="%m/%d/%y")
female_mortality$End.Age[female_mortality$Age==female_mortality$End.Age]<-female_mortality$End.Age[female_mortality$Age==female_mortality$End.Age]+0.1

models<-list()
dob.ranges$Uncertainty<-NA
for(i in 1:1000) {
  for(j in 1:nrow(dob.ranges)) {
    dob.ranges$Uncertainty[j]<-runif(n=1, min=dob.ranges$DownRange, max=dob.ranges$UpRange)
  }
  sim.data<-merge(female_mortality, dob.ranges, by="Individual")
  sim.data$Age<-round((sim.data$Age+sim.data$Uncertainty/365.25),1)
  sim.data$End.Age<-round((sim.data$End.Age+sim.data$Uncertainty/365.25),1)
  sim.data$End.Age[sim.data$Age==sim.data$End.Age]<-sim.data$End.Age[sim.data$Age==sim.data$End.Age]+0.1
  surv_object<-Surv(time=sim.data$Age, time2=sim.data$End.Age, event=sim.data$Censored)
  models[[i]]<-coxme(surv_object ~ poly(N.Females,2) + N.Males + (1|Unit), data = sim.data)
}

models
summary(model.avg(models))

# Extract effect sizes
effect.unit <- effects::effect(term= "poly(N.Females,2)", xlevels=10000, mod= mortality.model)
x_unit<-as.data.frame(effect.unit)

# Label extracted model fit for unit size
x_unit$unit.cat<-"SMALL"
x_unit$unit.cat[x_unit$N.Females>4.5]<-"MEDIUM"
x_unit$unit.cat[x_unit$N.Females>7.5]<-"LARGE"

## Make data points for Plot1A
unit.deaths<-female_mortality %>%
  group_by(round(N.Females)) %>%
  dplyr::summarise(N.Years=sum(End.Age-Age), N.Deaths=sum(Censored), Death.Rate=(N.Deaths/N.Years))
colnames(unit.deaths)<-c("N.Females", "N.Years", "N.Deaths", "Death.Rate")

## Make Plot 1A
library(scales)
plot1A<-ggplot(data = x_unit,aes(x=N.Females,y=fit)) + geom_ribbon(data= x_unit, aes(x=N.Females, ymin=lower, ymax=upper, fill=unit.cat), alpha= 0.8) +
  scale_fill_manual(values=c("#fec44f","#ec7014","#993404")) +
  geom_line(color="black") + theme(legend.position = "none") + ylim(0,.4) +
  geom_point(data=unit.deaths,aes(y=Death.Rate,x=N.Females, size=N.Years)) +
  scale_y_continuous(limits=c(0, 0.35), breaks=seq(0,0.3,0.1)) +
  xlab("# females") + ylab("deaths per female-year") + scale_x_continuous(limits=c(1,12), breaks=seq(2,12,2)) + theme(legend.position=c(0.8,0.8)) + labs(size="# female years") + theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)) +
  guides(size=guide_legend("# female years"), fill = FALSE) +
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1A

death.rate<-sum(female_mortality$Censored)/sum(female_mortality$End.Age-female_mortality$Age)
death.rate

## Make Plot 1B
mortality_data_summary<-ddply(female_mortality,. (unit.cat), summarize, ttl_fem_yrs=sum(End.Age-Age), ttl_deaths = sum(Censored), death_rate = ttl_deaths/ttl_fem_yrs, stdev = sqrt(death_rate*(abs(1-death_rate))/ttl_fem_yrs), avg_age=mean(Age), stdev_age=sd(Age))
mortality_data_summary$unit.cat<-factor(mortality_data_summary$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))
plot1B<-ggplot(mortality_data_summary, aes(x=as.factor(unit.cat),y=death_rate, fill=unit.cat)) + 
  geom_hline(yintercept=death.rate, color="#5a5a5a", linetype="dotted",size=0.8)+
  geom_errorbar(aes(ymin=death_rate-stdev, ymax=death_rate+stdev), width=.1, position=position_dodge(0.5)) +
  geom_point(shape=21,size=5, position=position_dodge(0.5), fill=c("#993404","#ec7014","#fec44f"), alpha=1) +
  ylab("")+
  scale_y_continuous(limits=c(0, 0.35), breaks=seq(0,0.3,0.1)) +
  xlab("")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1B

#############################################################################
###                     FEMALE REPRODUCTIVE PERFORMANCE                   ###                     
###                                 1-month                               ###
#############################################################################
# Construct model
repro.model<-glmer(Success ~ poly(N.Females,2) + poly(Age, 2) + N.Males + (1|Individual) + (1|Unit) + (1|Year), family="binomial", data=female_repro_1mos, control=glmerControl(optimizer="bobyqa")) 
summary(repro.model)
unique(female_repro_1mos$Unit)
cor.test(female_repro_1mos$N.Males, female_repro_1mos$N.Females, method="pearson")

## Extract model results
effect.unit <- effects::effect(term= "poly(N.Females,2)", xlevels=list(N.Females=seq(1,12,0.001)), mod= repro.model)
x_unit<-as.data.frame(effect.unit)
x_unit$unit.cat<-"SMALL"
x_unit$unit.cat[x_unit$N.Females>4.5]<-"MEDIUM"
x_unit$unit.cat[x_unit$N.Females>7.5]<-"LARGE"
female_repro_1mos$unit.cat<-"SMALL"
female_repro_1mos$unit.cat[female_repro_1mos$N.Females>4.5]<-"MEDIUM"
female_repro_1mos$unit.cat[female_repro_1mos$N.Females>7.5]<-"LARGE"
female_repro_1mos$unit.cat<-factor(female_repro_1mos$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))

## Build points for Plot 1C
infant.by.size<-female_repro_1mos %>%
  group_by(N.Females) %>%
  dplyr::summarise(Infants=sum(Success), N.Months=n(), N.Year=N.Months/12, Infants.Year=Infants/N.Year)

infant.by.size$unit.cat<-"SMALL"
infant.by.size$unit.cat[infant.by.size$N.Females>4.5]<-"MEDIUM"
infant.by.size$unit.cat[infant.by.size$N.Females>7.5]<-"LARGE"
infant.by.size$unit.cat<-factor(infant.by.size$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))

## Build histogram of unit sizes - Figure S1
infant.by.size$N.Females2<-as.character(infant.by.size$N.Females)
ggplot(data=infant.by.size, aes(x=N.Females, y=N.Year, fill=unit.cat)) +geom_bar(stat = "identity", width=1, color="black", alpha=0.8) +  scale_fill_manual(values=c("#993404","#ec7014","#fec44f")) +
  xlab("# females") + ylab("# female-years") + theme(legend.position = "none") + scale_x_continuous(breaks=seq(2,12,2), labels=seq(2,12,2)) + geom_vline(xintercept=mean(female_repro_1mos$N.Females), lty=2) +
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

# Make Plot 1C
plot1C<-ggplot(data = x_unit,aes(x=N.Females,y=fit*12)) + geom_ribbon(data= x_unit, aes(x=N.Females, ymin=lower*12, ymax=upper*12, fill=unit.cat), alpha= 0.8) +
  scale_fill_manual(values=c("#fec44f","#ec7014","#993404")) +
  geom_point(data=infant.by.size,aes(y=Infants.Year,x=N.Females, size=N.Year)) +
  xlab("# females") + ylab("surviving offspring\n per female-year") + scale_x_continuous(limits=c(1,12), breaks=seq(2,12,2))+ geom_line() + theme(legend.position=c(0.85,0.85)) + labs(size="# female-years") + theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)) +
  guides(size=guide_legend("# female years"), fill = FALSE) +
  scale_y_continuous(limits=c(0,0.42),breaks=seq(0,0.40,by=0.10)) + theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1C

## Calculate success rate
success_rate<-sum(female_repro_1mos$Success)/(nrow(female_repro_1mos)/12)
success_rate

## Make Plot 1D
repro_data_summary<-ddply(female_repro_1mos,. (unit.cat), summarize, ttl_fem_yrs=length(Individual), ttl_success = sum(Success), success_rate = ttl_success/ttl_fem_yrs, stdev = sqrt(success_rate*(abs(1-success_rate))/ttl_fem_yrs), avg_age=mean(Age), stdev_age=sd(Age))
plot1D<-ggplot(repro_data_summary, aes(x=as.factor(unit.cat),y=success_rate*12, fill=unit.cat)) + 
  geom_hline(yintercept=success_rate, color="#5a5a5a", linetype="dotted",size=0.8)+
  geom_errorbar(aes(ymin=success_rate*12-stdev*12, ymax=success_rate*12+stdev*12), width=.1, position=position_dodge(0.5)) +
  geom_point(shape=21,size=5, position=position_dodge(0.5), fill=c("#993404","#ec7014","#fec44f"), alpha=1) +
  ylab("") +
  scale_y_continuous(limits=c(0,0.42),breaks=seq(0,0.40,by=0.10)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1D

## CALCULATE RS AT EACH CATEGORY
(16.0-6.06)*(0.02261307)*(12)
(19.7-6.06)*(0.02658396)*(12)
(16.0-6.06)*(0.01901408)*(12)

plot1A+plot1B+plot1C+plot1D

#############################################################################
###                     FEMALE REPRODUCTIVE PERFORMANCE                   ###                     
###                                  6-month                              ###
#############################################################################
# Construct model
repro.model.6mo<-glmer(Success ~ poly(N.Females,2) + poly(Avg.Age,2) + N.Males + offset(N.Month) + (1|Individual) + (1|Unit) + (1|SuccessBin), family="binomial", data=female_repro_6mos, control=glmerControl(optimizer="bobyqa")) 
summary(repro.model.6mo)

## Extract model results
effect.unit <- effects::effect(term= "poly(N.Females,2)", xlevels=100, mod= repro.model.6mo)
x_unit<-as.data.frame(effect.unit)
ggplot(data = x_unit,aes(x=N.Females,y=fit*12)) + geom_line() + geom_ribbon(data= x_unit, aes(x=N.Females, ymin=lower*12, ymax=upper*12), alpha= 0.3, fill="blue") +
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

female_repro_6mos$unit.cat<-"SMALL"
female_repro_6mos$unit.cat[female_repro_6mos$N.Females>4]<-"MEDIUM"
female_repro_6mos$unit.cat[female_repro_6mos$N.Females>7]<-"LARGE"
female_repro_6mos$unit.cat<-factor(female_repro_6mos$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))
rate<-sum(female_repro_6mos$Success)/nrow(female_repro_6mos)

## Construct plot comparable to Plot 1D
mortality_data_summary<-ddply(female_repro_6mos,. (unit.cat), summarize, ttl_fem = mean(N.Females), ttl_fem_yrs=sum(N.Month)/6, ttl_success = sum(Success), success_rate = ttl_success/ttl_fem_yrs, stdev = sqrt(success_rate*(abs(1-success_rate))/ttl_fem_yrs), avg_age=mean(Avg.Age), stdev_age=sd(Avg.Age))
plot.6mo<-ggplot(mortality_data_summary, aes(x=as.factor(unit.cat),y=success_rate*2, fill=unit.cat)) + 
  geom_hline(yintercept=rate*2, color="#5a5a5a", linetype="dotted",size=0.8)+
  geom_errorbar(aes(ymin=success_rate*2-stdev, ymax=success_rate*2+stdev), width=.1, position=position_dodge(0.5)) +
  geom_point(shape=21,size=5, position=position_dodge(0.5), fill=c("#993404","#ec7014","#fec44f"), alpha=1) +
  ylab("\nproportion of 6-month periods successful")+
  scale_y_continuous(limits=c(0,0.4),breaks=seq(0,0.4,by=0.1)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot.6mo

#############################################################################
###                         INTERBIRTH INTERVALS                          ###
#############################################################################
ibi_data$unit.cat<-"SMALL"
ibi_data$unit.cat[ibi_data$N.Females>4]<-"MEDIUM"
ibi_data$unit.cat[ibi_data$N.Females>7]<-"LARGE"

# Construct model with all data
IBI.model<-lmer(data=ibi_data, IBI ~ poly(N.Females,2) + poly(Mom.Age,2) + Sex + (1|Individual) + (1|Unit), control=lmerControl(optimizer="bobyqa"))
summary(IBI.model)

## GET RID OF TAKEOVERS
ibi_no_TO<-ibi_data[ibi_data$Takeover=="N",]

# Construct model with takeovers removed
model.NO<-lmer(data=ibi_no_TO, IBI ~ poly(N.Females,2) + poly(Mom.Age,2) + Sex + (1|Individual) + (1|Unit), control=lmerControl(optimizer="bobyqa"))
summary(model.NO)

## Plot both results as beeswarms
theme_set(theme_minimal(base_size=12))

std <- function(x) sd(x)/sqrt(length(x))

ibi_cat1<-ibi_data %>%
  group_by(unit.cat) %>%
  dplyr::summarise(Mean.IBI=mean(IBI), min=Mean.IBI-1.96*std(IBI), max=Mean.IBI+1.96*std(IBI), N.IBI=n())

## Plot 2A and 2B
## Extract model results
effect.unit <- effects::effect(term= "poly(N.Females,2)", xlevels=list(N.Females=seq(1,10,0.001)), mod= IBI.model)
x_unit<-as.data.frame(effect.unit)
x_unit
x_unit$unit.cat<-"SMALL"
x_unit$unit.cat[x_unit$N.Females>=4.5]<-"MEDIUM"
x_unit$unit.cat[x_unit$N.Females>=7.5]<-"LARGE"

plot2A<-ggplot() + scale_color_manual(values=c("#993404","#ec7014","#fec44f"))+
  ylab("interbirth interval (years)") + xlab("") + ylim(0.5,5.5) + scale_x_continuous(limits=c(1,10),breaks=seq(2,10,2)) +
  geom_ribbon(data= x_unit, aes(x=N.Females, ymin=lower/365.25, ymax=upper/365.25, fill=unit.cat), alpha= 0.8) +
  geom_line(data=x_unit, aes(x=N.Females, y=fit/365.25)) + theme_minimal(base_size = 12) +
  geom_point(data = ibi_data, aes(x=N.Females, y=IBI/365.25),size=1.3) + theme(legend.position="none") +
  scale_fill_manual(values=c("#fec44f","#ec7014","#993404")) +
  annotate("text", x = 5.5, y = 5.2, label = "P=0.043", fontface=2) +
  labs(title="with takeovers") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot2A

plot2B<-ggplot() + scale_color_manual(values=c("#993404","#ec7014","#fec44f")) + geom_point(data = ibi_no_TO, aes(x=N.Females, y=IBI/365.25),size=1.3) + theme(legend.position="none") +
  ylab("interbirth interval (years)") + xlab("") + ylim(0.5,5.5) + scale_x_continuous(limits=c(1,10),breaks=seq(2,10,2)) +
  annotate("text", x = 5.5, y = 5.2, label = "NS", fontface=2) + xlab("# females") + theme_minimal(base_size = 12) +
  labs(title="without takeovers") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot2B


#############################################################################
###                            INFANT SURVIVAL                            ###
#############################################################################
theme_set(theme_minimal(base_size=20))

surv_object<-Surv(time=infant_survival$End.Age, event=infant_survival$Censored)
fit.null <- survfit(surv_object ~ 1, data = infant_survival)
fit.unit.cat <- survfit(surv_object ~ unit.cat, data = infant_survival)

# Construct model with all data
coxme.infant <- coxme(surv_object ~ poly(N.Females,2) + N.Males + Sex + (1|Mom) + (1|Unit) + (1|Year),  data = infant_survival)
summary(coxme.infant)

length(unique(infant_survival$Unit))
length(unique(infant_survival$Mom))
length(unique(infant_survival$Year))

# Plot 2C - survival curve
infant_survival$unit.cat<-factor(infant_survival$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))

sum(infant_survival$Censored[infant_survival$unit.cat=="SMALL"])/nrow(infant_survival[infant_survival$unit.cat=="SMALL",])
sum(infant_survival$Censored[infant_survival$unit.cat=="LARGE"])/nrow(infant_survival[infant_survival$unit.cat=="LARGE",])

library(RColorBrewer)
library(survminer)
palette1<-brewer.pal(3, name = "Dark2")
survival2<-ggsurvplot(fit.unit.cat, data = infant_survival, xlim=c(0,1.5), ylim=c(.5,1), legend="none", conf.int = F,  censor=F, xlab="age (years)", ylab="proportion infants surviving", legend.labs=c("LARGE", "MEDIUM", "SMALL"), legend.title="", break.time.by=0.5, palette = c("#fec44f","#ec7014","#993404"), font.x=c(size=20), font.y=c(size=20))
survival2
plot2C<-survival2$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot2C 

# Plot null survival

survival.null<-ggsurvplot(fit.null, data = infant_survival, xlim=c(0,1.5), ylim=c(0,1), legend="none",  censor=F, xlab="Age (years)", ylab="Proportion infants surviving", legend.title="", break.time.by=0.5, palette = palette1)
null.plot<-survival.null$plot +  theme(plot.title = element_text(hjust = 0,face = "bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
null.plot

# Combine plots for Figure 2 panels

layout <- "
ACC
BCC
"

plot2A
plot2B
plot2C 

library(patchwork)
plot2A + plot2B + plot2C + plot_layout(design=layout) + plot_annotation(tag_levels = "a") 

#############################################################################
###                             TAKEOVER RATES                            ###
#############################################################################
# Construct model
unit_takeovers$N.Year<-unit_takeovers$N.Months/12
unit_takeovers$Male.Cat<-"No Followers"
unit_takeovers$Male.Cat[unit_takeovers$N.Males>1]<-"Followers"
mod4<-glmer(Takeovers ~ poly(N.Females,2) + N.Males + offset(N.Year) + (1|Unit)+(1|Year), family="poisson", data=unit_takeovers, control=glmerControl(optimizer="bobyqa")) 
summary(mod4)

## Calculate takeover rates
unit_takeovers$unit.cat<-factor(unit_takeovers$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))
takeover_rate<-sum(unit_takeovers$Takeovers)/sum(unit_takeovers$N.Months/12)
1/(sum(unit_takeovers$Takeovers[unit_takeovers$unit.cat=="SMALL"])/sum(unit_takeovers$N.Months[unit_takeovers$unit.cat=="SMALL"]))/12
1/(sum(unit_takeovers$Takeovers[unit_takeovers$unit.cat=="MEDIUM"])/sum(unit_takeovers$N.Months[unit_takeovers$unit.cat=="MEDIUM"]))/12
1/(sum(unit_takeovers$Takeovers[unit_takeovers$unit.cat=="LARGE"])/sum(unit_takeovers$N.Months[unit_takeovers$unit.cat=="LARGE"]))/12

# PLOT 3A
takeover_data_summary<-ddply(unit_takeovers,. (unit.cat), summarize, ttl_fem_yrs=sum(N.Months)/12, ttl_takeovers = sum(Takeovers), takeover_rate = ttl_takeovers/ttl_fem_yrs, stdev = sqrt(takeover_rate*(abs(1-takeover_rate))/ttl_fem_yrs))
plot3A<-ggplot(takeover_data_summary, aes(x=as.factor(unit.cat),y=takeover_rate)) +
  geom_hline(yintercept=takeover_rate, color="#5a5a5a", linetype="dotted",size=0.8)+ 
  geom_errorbar(aes(ymin=takeover_rate-stdev, ymax=takeover_rate+stdev), width=.1, position=position_dodge(0.5)) + 
  geom_point(size=5, position=position_dodge(0.5), shape=21, fill=c("#993404", "#ec7014","#fec44f"), alpha=1) +
  ylab("takeovers per unit-year")+ theme(legend.position = "none") +
  scale_y_continuous(limits=c(0,.75),breaks=seq(0,1.25,by=0.25),labels=c(0,0.25,0.50,0.75,1, 1.25)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot3A

unit_takeovers$male.cat<-"NO FOLLOWERS"
unit_takeovers$male.cat[unit_takeovers$N.Males>1]<-"FOLLOWERS"
unit_takeovers$male.cat<-factor(unit_takeovers$male.cat, levels=c("NO FOLLOWERS", "FOLLOWERS"))


#############################################################################
###                        CAUSES OF INFANT MORTALITY                     ###
#############################################################################
## Model 1: Infanticide
model.infanticides<-glmer(data=infant_mortality_causes, Infanticide ~ poly(N.Females,2) + (1|Unit), family="binomial")
summary(model.infanticides)

# Model 2: Maternal death
model.momdeath<-glmer(data=infant_mortality_causes, MaternalDeath ~ poly(N.Females,2) + (1|Unit), family="binomial")
summary(model.momdeath)

# Model 3: Unknown deaths
model.unknown<-glmer(data=infant_mortality_causes, Unknown ~ poly(N.Females,2) + (1|Unit), family="binomial")
summary(model.unknown)

infant.cause.summary<-infant_mortality_causes %>%
  group_by(unit.cat) %>%
  dplyr::summarise(Infanticide=sum(Infanticide)/n(), MatDeath=sum(MaternalDeath)/n(), Unknown=sum(Unknown)/n())

infant.cause.summary <- melt(infant.cause.summary, id=c("unit.cat"))
infant.cause.summary$unit.cat<-factor(infant.cause.summary$unit.cat, levels=c("SMALL", "MEDIUM", "LARGE"))

# Plot 3B
palette1<-brewer.pal(name="Dark2",8)
palette1<-palette1[c(1,3,8)]
plot3B<-ggplot(data=infant.cause.summary,aes(x=unit.cat, y=value, fill=variable)) + 
  geom_bar(color="black",stat="identity") +
  xlab("unit size category") + ylab("infant mortality")  + scale_fill_manual(name = "Likely cause of\ninfant death", labels=c("Infanticide", "Maternal death", "Unknown"),values=palette1) +
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot3B

