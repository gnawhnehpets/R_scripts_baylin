rm(list=ls())
for(i in 1:40){gc()}
library(KMsurv)
library(survival)
#setwd("Y:./././users/shwang26/m084b trial/")
setwd("/home/steve/.gvfs//onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/")
#pre <- read.table("patientstatus_pre_allsamples.csv", header = TRUE)
both <- read.table("patientstatus_prepost_allsamples.csv", header=TRUE)
both <- read.table("patientstatus_genes_allsamples.csv", header=TRUE) #same as prepost_allsamples but with new genes
both
# both <- read.table("patientstatus_genes_pairsonly.csv", header=TRUE) #pairs only
# post <- read.table("patientstatus_post.csv")
both
dev.off()
pre <- both[1:34,]
pre
post <- both[35:51,]
post
###############################
table <- pre #pre, post, both
table
###############################
pfsdays <- table[,4]
osdays <- table[,6]
posttrtdays <- osdays-pfsdays
statusearly <- table[,8]
statustwomo <- table[,9]
statusfivesix <- table[,10]
statusoverall <- table[,11]
statusbest <- table[,12]
ecad <- table[,13]
p16 <- table[,14]
tms1 <- table[,15]
timp1 <- table[,16]
timp2 <- table[,17]
timp3 <- table[,18]
timp4 <- table[,19]
mlh1 <- table[,20]
mgmt <- table[,21]
pd1 <- table[,22]
pdl1 <- table[,23]
pdl2 <- table[,24]
sfrp1 <- table[,25]
tfp1 <- table[,26]
gata4 <- table[,27]
gata5 <- table[,28]
apc <- table[,29]
chfr <- table[,30]
rassf1 <- table[,31]
hin1 <- table[,32]
group <- table[,33] #1: pre, #2: post

########################
# CHANGE DESIGNATION PER GROUP
days2 <- osdays
status2 <- statusoverall
############################################################################################################################################################################################################################################################
# COX PROPORTIONAL HAZARDS
############################################################################################################################################################################################################################################################
# 2-group comparison (http://www.ics.uci.edu/~vqnguyen/stat255/Lecture06.pdf)
# Question: how does (gene) expression (1: low expression, 2: high expression) effect survival?
fit <- coxph(Surv(days2, status2) ~ ecad, data = table)
fit2 <- coxph(Surv(days2, status2) ~ p16, data = table)
fit3 <- coxph(Surv(days2, status2) ~ tms1, data = table)
fit4 <- coxph(Surv(days2, status2) ~ timp1, data = table)
fit5 <- coxph(Surv(days2, status2) ~ timp2, data = table)
fit6 <- coxph(Surv(days2, status2) ~ timp3, data = table)
fit7 <- coxph(Surv(days2, status2) ~ timp4, data = table)
fit8 <- coxph(Surv(days2, status2) ~ mlh1, data = table)
fit9 <- coxph(Surv(days2, status2) ~ mgmt, data = table)
fit10 <- coxph(Surv(days2, status2) ~ pd1, data = table)
fit11 <- coxph(Surv(days2, status2) ~ pdl1, data = table)
fit12 <- coxph(Surv(days2, status2) ~ pdl2, data = table)
fit13 <- coxph(Surv(days2, status2) ~ sfrp1, data = table)
fit14 <- coxph(Surv(days2, status2) ~ tfp1, data = table)
fit15 <- coxph(Surv(days2, status2) ~ gata4, data = table)
fit16 <- coxph(Surv(days2, status2) ~ gata5, data = table)
fit17 <- coxph(Surv(days2, status2) ~ apc, data = table)
fit18 <- coxph(Surv(days2, status2) ~ chfr, data = table)
fit19 <- coxph(Surv(days2, status2) ~ rassf1, data = table)
fit20 <- coxph(Surv(days2, status2) ~ hin1, data = table)
exp(fit$coefficient)

# Conclusion: The risk of death is exp(fit$coefficient) times higher for the (2) as compared to the (1) #table = both

summary(fit)
#exp(coef): 1.74/2.69
#z: 1.53/1.6
#p: .126/.107
#95 CI: contains 1  x2
summary(fit2)
#exp(coef): .77/.4565
#z: -.74/-1.4
#p: .455/.162
#95 CI: contains 1  x2
summary(fit3)
#exp(coef): .87/.543
#z: -.37/-1.18
#p: .71/.24
#95 CI: contains 1  x2
summary(fit4) #
#exp(coef): 1.23/.3582
#z: .60/-1.8
#p: .55/.0735
#95 CI: contains 1 x2
summary(fit5)
#exp(coef): .83/.804
#z: -.5/-.42
#p: .61/.68
#95 CI: contains 1 x2
summary(fit6)
#exp(coef): .88/2.16
#z: -.36/1.44
#p: .72/.147
#95 CI: contains 1 x2
summary(fit7) #****pre/timp4
# Conclusion: The risk of death is 2.17 times higher for the pre-samples with high timp4 expression as compared to the pre-samples with low timp4 expression.
#exp(coef): 2.17/1.4
#z: 2.047/ .65
#p: .0407/ .516
#95 CI: does not contain 1 / contains 1
summary(fit8)
#exp(coef): 1.55/.91
#z: 1.2/-.172
#p: .233/.86
#95 CI: contains 1  x2
summary(fit9)
#exp(coef): 1.11/.5
#z: .301/-1.3
#p: .763/.198
#95 CI: contains 1  x2
summary(fit10)
#exp(coef): .89/.66
#z: -.34/-.8
#p: .7/.43
#95 CI: contains 1  x2
summary(fit11)
#exp(coef): 1.1/1.28
#z: .18/.48
#p: .85/.63
#95 CI: contains 1  x2
summary(fit12)
#exp(coef): 1.0/1.44
#z: .003/.66
#p: .998/.51
#95 CI: contains 1  x2
summary(fit13)
#exp(coef): .75
#z:  -.77
#p: .44
#95 CI:  contains 1
summary(fit14)
#exp(coef): .61
#z: -1.3
#p: .18
#95 CI: contains 1
summary(fit15)
#exp(coef): .775
#z:  -.73
#p: .47
#95 CI: contains 1 
summary(fit16)
#exp(coef): .82
#z: -.58
#p: .56
#95 CI:  contains 1
summary(fit17)
#exp(coef): .84
#z: -.47
#p: .64
#95 CI: contains 1
summary(fit18)
#exp(coef): 1.1
#z: .09
#p: .93
#95 CI: contains 1
summary(fit19)
#exp(coef): .81
#z: -.59
#p: .56
#95 CI: contains 1 
summary(fit20)
#exp(coef): 1.2
#z: .512
#p: .61
#95 CI: contains 1 


# > summary(fit)
# Call:
#   coxph(formula = Surv(days2, status2) ~ ecad, data = table)
# 
# n= 51, number of events= 51 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# ecad 0.6481    1.9120   0.3030 2.139   0.0324 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# ecad     1.912      0.523     1.056     3.463
# 
# Concordance= 0.572  (se = 0.042 )
# Rsquare= 0.087   (max possible= 0.997 )
# Likelihood ratio test= 4.62  on 1 df,   p=0.03168
# Wald test            = 4.58  on 1 df,   p=0.03243
# Score (logrank) test = 4.71  on 1 df,   p=0.02999

# Conclusion: The risk of death is exp(fit$coefficient) times higher for the (1) as compared to the (2) #table = both
# z = 2.14, p = .0324
# hypothesis Ho = beta = 0 (exp(beta)=1)
# 95% CI for beta= [1.056, 3.463] #does not contain 1 --> effect low/high expression of (gene) is significantly different than 1 @5% level

############################################################################################################################################################################################################################################################
# Multiple regression #stratified levels/multiple groups/stages of disease 
# pre-group vs. post-group 
# beta1: log-relative hazard (hazard ratio) comparing two groups that differ in treatment received. 
# exp(beta1): hazard ratio comparing two groups that differ in treatment received 
# beta1: effect of treatment received adjusting 
############################################################################################################################################################################################################################################################
# multiple regression
# What is the interpretation of exp(coef)? hazard ratio comparing a group that showed stable disease (SD) during the trial to a group that never showed stable disease (PD) where both groups have high or low expression of ecad
# What is the interpretation of exp(coef)? hazard ratio comparing the pre/post group, where both groups showed high or low expression of (gene)
# Does not matter what ecad expression the two groups have (just that they be the same) - the effect of being in pre or post group is the same
# Model assumes the effect of AZA is the same regardless of (gene) expression level
# Model assumes the effect of (gene) expressional level is the same regardless of being treated with AZA (post) or not (pre)
#fit <- coxph(Surv(time, delta) ~ age + factor(stage), data=larynx)

gene <- pdl2 #ecad, p16, tms1, timp1, timp2, timp3, timp4, mlh1, mgmt, pd1, pdl1, pdl2, sfrp1, tfp1, gata4, gata5, apc, chfr, rassf1, hin1
gname <- "pdl2"
#cox <- coxph( Surv( days2, status2 ) ~ gene + factor(statusbest), data=table)
cox <- coxph( Surv( days2, status2 ) ~ gene + factor(group), data=table) # if using 'both'
cox <- coxph( Surv( days2, status2 ) ~ gene, data=table) # if using 'pre'
summary(cox)

# Interpretation
# http://www.ics.uci.edu/~vqnguyen/stat255/Lecture06.pdf
# ecad:  Estimate the risk of death to be 1.9x higher in post samples than that of pre samples with respect to ecad expression # p = .032  
# timp4: Estimate the risk of death to be 1.96x higher in post samples than that of pre samples with respect to timp4 expression # p = .024
# p = .032  

title <- paste("overall survival between patients based on", gname, "expression", sep=" ")
xtitle <- "Time from study start (days)"
f <- survfit( Surv( days2, status2 ) ~ gene, data=table)
cox2 <- survfit(cox)

# Alright so it seems like majority has voted for Option B. Tony, you're planning to meet us there, right? Michael, any update on your decision?
# Here's the tentative schedule.
# Thursday: I fly in 930p Thursday night, and I plan to drink some beer, go through the gear checklist, and prep food for the trip. Everyone is more than welcome to crash at my parents' house for the night.
# Friday: Campsites are first-come, first-served for Friday night so it would be best to leave in the morning. Drive is ~6 hours, but we'll stop by Louisville for lunch. Let's aim to arrive before 3pm.
# Saturday: Canoe company has us booked for a 10am trip, which means we gotta pack up @9a. We'll pull off on Sand Cave Island for lunch, canoe/fish down the river, then set up camp @ Crump Island for the night. 
# Sunday: Wake up, pack up, canoe the rest of the route, get shuttled back to campground, cave tour, then head home.
# 
# Here's a Google map of our itinerary: https://mapsengine.google.com/map/edit?mid=zoqBJ8pwd3GQ.kM4ONQ9pzqg0

# Please be familiar with the 'Bearmuda' principle of setting up camp <http://tinyurl.com/onbezr7>. Bear canisters are not required at this national park, but bears are still a potential risk. 

# Please become familiar with the rules of the park: http://www.nps.gov/maca/planyourvisit/loader.cfm?csModule=security/getfile&PageID=107668
 


plot(f, col=1:2, lty=1:2, main = title, xlab=xtitle, ylab="Survival probability", conf.int=TRUE, xlim = c(0,400))
polygon(c(mod$time, rev(mod$time)), c(mod$upper, rev(mod$lower)), col = 'lightgrey', border = NA)
lines(mod$time, mod$upper, type='l', lty=2)
lines(mod$time, mod$lower, type='l', lty=2)
lines(mod$time, mod$surv, type="l", lty=1)
mod <- summary(f)
class(mod)

legend( "topright", col=1:2, lty=1:2, legend=c("low exp", "high exp"), bty="n" )
dev.copy(png, width=1000, height=1000, paste("os_exp_nocox_",gname,"-pre.png",sep=""))
dev.off()

plot(cox2, col=1:2, lty=1:2, main = title, xlab=xtitle, ylab="Survival probability", conf.int=FALSE, xlim = c(0,400))
legend( "topright", col=1:2, lty=1:2, legend=c("low exp", "high exp"), bty="n" )
dev.copy(png, width=1000, height=1000, paste("os_exp_coxph_",gname,"-pre.png",sep=""))
dev.off()
#for(i in 1:10){dev.off()}

#abline(h=.5, col=222, lty=3)
############################################################################################################################################################################################################################################################
# log-rank test
# ecad diff
x <- statusbest
#ecad, p16, tms1, timp1, timp2, timp3, timp4, mlh1, mgmt, pd1, pdl1, pdl2
survdiff(Surv(days2, status2) ~ x, data=table)
# timp4: survival times different between patients with low timp4 expression and high timp4 expression #pre

############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################
#setwd("Y:./././users/shwang26/m084b trial/")
setwd("/home/steve/.gvfs//onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/m084b trial/")
meta <- read.table("patientstatus_statusbest.csv", header=TRUE)
meta
#remove 1644 outlier
table <- meta[-47,]
# table <- meta
table

pfsdays <- table[,4]
osdays <- table[,6]
posttrtdays <- osdays-pfsdays
statusearly <- table[,8]
statustwomo <- table[,9]
statusfivesix <- table[,10]
statusoverall <- table[,11]
statusbest <- table[,12]
statusbest
#######################

title <- "survival of patients with colon cancer (osdays/statusoverall)"
xtitle <- "time since start of trial (in days)"
days2 <- osdays
status2 <- statusoverall

title <- "progression-free survival in trial (pfsdays/statusoverall)"
xtitle <- "time since start of trial (in days)"
days2 <- pfsdays
status2 <- statusoverall

# This wouldn't make sense since all of them deceased
# title <- "progression-free survival in trial (pfsdays/statusfivesix)"
# xtitle <- "time since start of trial (in days)"
# days2 <- pfsdays
# status2 <- statusfivesix
# 
# title <- "progression-free survival in trial (pfsdays/statustwo)"
# xtitle <- "time since start of trial (in days)"
# days2 <- pfsdays
# status2 <- statustwomo

title <- "survival of patients post treatment (posttrtdays/statusoverall)"
xtitle <- "time since end of trial (in days)"
days2 <- posttrtdays
status2 <- statusoverall

#######################
mydat
mydat <- data.frame(cbind(days2, status2))
mydat
mysurv <- with(mydat, Surv(days2, status2))
mysurv
f <- survfit(mysurv ~ 1,type="kaplan-meier",data=mydat)
f[2]
names(f)
f$n.event

plot(f, conf.int=TRUE, main = title, xlab= xtitle)
abline(v=60, col=20, lty=3)
plot(f, conf.int=FALSE, main = title, xlab= xtitle)
abline(v=55, col=20, lty=3)

abline(v=c(60,120,180,240,300,360,420,480,540,600,660,720,780,840,900,960,1020),lty=3)
abline(h=.5, col=222, lty=3)

########################
# compare with 4 patients who did not partcipate in trial
days2 <- osdays
status2 <- statusoverall
status2 <- statusbest

# plot(survfit( Surv( time, death ) ~ im, data=btrial), lty=1:2
#      , xlab="Time from study start (months)"
#      , ylab="Survival probability")
# legend(  80, 1, lty=1:2, legend=c("first", "second"), bty="n" )
#im <- table[,8]
im <- table[,12]
im
mydat <- data.frame(cbind(days2, status2, im))
mydat
mysurv <- with(mydat, Surv(days2, status2) ~ im)
mysurv
f <- survfit(mysurv,type="kaplan-meier",data=mydat)
f[1]
f[2]
title <- "overall survival between patients based on best response"
xtitle <- "Time from study start (days)"
f <- survfit( Surv( days2, status2 ) ~ im, data=mydat)
plot(f, col=1:3, lty=1:3, main = title, 
     xlab=xtitle, ylab="Survival probability", conf.int=FALSE, xlim = c(0,400))
legend( "topright", col=1:3, lty=1:3, legend=c("NE/NA", "PD", "SD"), bty="n" )
abline(h=.5, col=222, lty=3)

# summary(f)
# mod <- summary(f)
# with(mod,plot(day2,surv,type="n"))
# with(mod,polygon(c(day2,rev(time)),c(lower,rev(upper)),
#                  col = "grey75", border = FALSE))
# with(mod,lines(day2),surv,type="s"))


x <- c(0:5)
y.a <- seq(from=1,to=.5,length=6)
y.b <- c(1,.7,.56,.52,.52,.5)

plot(c(0,5),c(0,1),type='n',xlab='Follow-up Time (years)',ylab='Survival probability',main='Survival experience for groups A and B')
lines(x,y.a,col='red',lwd=2)
points(x,y.a,col='red',cex=1.5)
lines(x,y.b,col='blue',lwd=2)
points(x,y.b,col='blue',cex=1.5)

dev.off()
plot(c(0,12),c(1,6),axes=F,xlab="Weeks",ylab="Patients",type="n")
axis(side=2,at=c(1:6),labels=c("P1","P2","P3","P4","P5","P6"))
axis(side=1,at=seq(0,12,2),labels=c("0","2","4","6","8","10","12"))

lines(c(0:5),rep(6,6)); points(5,6,cex=1.5,pch="X");
#P5
lines(c(0:12),rep(5,13)); points(12,5,cex=1.5,pch="O");
#P4
lines(c(4:6),rep(4,3)); points(6,4,cex=1.5,pch="O");
#P3
lines(c(5:12),rep(3,8)); points(12,3,cex=1.5,pch="O");
#P2
lines(c(4:10),rep(2,7)); points(10,2,cex=1.5,pch="O");
#P1
lines(c(9:12),rep(1,4)); points(12,1,cex=1.5,pch="X");
abline(v=12,lty=2)

dev.off()
plot(c(0,8),c(1,38),axes=F,xlab="Weeks",ylab="Patients",type="n")
patient <- NULL
for(i in 1:38){
  patient = c(patient, paste("P",i,sep=""))
}
patient

axis(side=2,at=c(1:38),labels=patient)
weeks <- as.character(c(1:8))
weeks
axis(side=1,at=seq(0,8,2),labels=weeks)
