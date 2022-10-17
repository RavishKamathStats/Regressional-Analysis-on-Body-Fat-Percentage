install.packages("car")
install.packages("tidyverse")
install.packages("matlib")
install.packages("leaps")
install.packages("MASS")
install.packages("ggplot2")
install.packages("qqplotr")
install.packages("SignifReg")
install.packages("moments")
install.packages("corrplot")
install.packages("coefplot")
install.packages("olsrr")
install.packages("qpcR")
library(tidyverse)
library(car)
library(matlib)
library(leaps)
library(MASS)
library(ggplot2)
library(qqplotr)
library(SignifReg)
library(corrplot)
library(coefplot)
library(moments)
library(olsrr)
library(qpcR)
attach(df0) 
plot(df0)

sink("./mod5.txt", append =  T)
sink()


#FUNCTIONS##################################################################
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}
model_fit_stats <- function(linear.model) {
  r.sqr <- summary(linear.model)$r.squared
  adj.r.sqr <- summary(linear.model)$adj.r.squared
  pre.r.sqr <- pred_r_squared(linear.model)
  PRESS <- PRESS(linear.model)
  return.df <- data.frame(r.squared = r.sqr, adj.r.squared = adj.r.sqr, pred.r.squared = pre.r.sqr, press = PRESS)
  return(return.df)
}

#DATA TABLES##################################################################
#Remove certain observations that are impossible:
view(df0)
df0 %>% filter(BodyFat <= 2)
df0 %>% filter(Height == 29.5)
df1 = df0 %>% filter(Height != 29.5) %>% filter(BodyFat >=2)
view(df1)

#Let us categorize Age and combine Height and Weight into BMI
BMI = (df1$Weight/(df1$Height)^2)*703
BMI = round(BMI, digits = 2)
catAge = cut(df1$Age, breaks = c(20,40,60,90), labels = c("Young", "Middle", "Older"), right = FALSE)
levels(catAge)
df2 = cbind(df1, BMI, catAge)
df2 = subset(df2, select= -c(Weight, Height, Age))
view(df2)


#Influential point removal:
df2 %>% filter(df2$BMI > 48)
df3_ob39rem = df2 %>% filter(BMI <= 41)
view(df3_ob39rem)

#FULL MODEL##################################################################
#Full MODEL
summary(lm(BodyFat ~ Thigh, data = df2))

fit0 = lm(BodyFat ~ ., data = df2)
summary(fit0)
anova(fit0)
residualPlots(fit0, type = "rstudent")
par(mfrow = c(2,1))
hist(fit0$resid, main = "Normalilty Plot for Body Fat Percentage", xlab = "Fitted Values")
t_i = rstudent(fit0) 
qqnorm(t_i)                                         
qqline(t_i, col = 2, lwd = 2)
skewness(fit0$residuals)
summary(influence.measures(fit0))
head(fortify(fit0))
vif(fit0)

residPlot = ggplot(aes(x = .fitted, y = .resid), data = fit0) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot

#Full model without observation 39
alt_fit0 = lm(BodyFat ~ ., data = df3_ob39rem)
summary(alt_fit0)
anova(alt_fit0)
residualPlots(alt_fit0, type = "rstudent")
par(mfrow = c(2,1))
hist(alt_fit0$resid, main = "Normalilty Plot for Body Fat Percentage", xlab = "Fitted Values")
alt_t_i = rstudent(alt_fit0) 
qqnorm(alt_t_i)                                         
qqline(alt_t_i, col = 2, lwd = 2)
skewness(alt_fit0$residuals)
summary(influence.measures(alt_fit0))
head(fortify(alt_fit0))
residPlot = ggplot(aes(x = .fitted, y = .resid), data = alt_fit0) 
+ geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
vif(alt_Fit0)

residPlot = ggplot(aes(x = .fitted, y = .resid), data = alt_fit0) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot



#MODEL BUILDING###########################################################

#Forward Selection based on p value
forward.fit <- ols_step_forward_p(fit0, penter = 0.25, details = TRUE)
fwdModFit = lm(BodyFat ~ Abdomen + Wrist + catAge + Neck + Hip + Forearm + Thigh, data = df2)

forward.altFit = ols_step_forward_p(alt_fit0, penter = 0.25, details = TRUE)
fwdModAltFit = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest + Hip + Neck + Forearm + Thigh, data = df3_ob39rem)
summary(fwdModAltFit)
anova(fwdModAltFit)

#Backward selection based on p-value
backward.fit <- ols_step_backward_p(fit0, prem = 0.1, progress = T , details = T )
bwdModFit = lm(BodyFat ~ Neck + Abdomen + Hip + Thigh + Forearm + Wrist + catAge, data = df2)

backward.altFit = ols_step_backward_p(alt_fit0, prem = 0.1, progress = T , details = T )
bwdModAltFit = lm(BodyFat ~ Chest + Abdomen + Wrist + BMI + catAge, data = df3_ob39rem )
summary(bwdModAltFit)
anova(bwdModAltFit)

#Stepwise Regression based on P-value
step.fit = ols_step_both_p(fit0, pent = 0.25, prem = 0.1, progress = TRUE , details = TRUE)
stepModFit = lm(BodyFat ~ Abdomen + Wrist + Hip + catAge + Forearm + Neck + Thigh, data = df2)

step.altFit = ols_step_both_p(alt_fit0, pent = 0.25, prem = 0.1, progress = TRUE , details = TRUE)
stepModAltFit = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest, data = df3_ob39rem)
summary(stepModAltFit)
anova(stepModAltFit)

# Stepwise, finding the best subset of the full Model
sub.altFit = ols_step_best_subset(alt_fit0)
sub.altFit

#MODELS TO CHOOSE###########################################################
# REDUCED MODELS 

#1st Model Regressors are Abdomen, catAge, BMI, Neck, Thigh, Hip, Forearm, Wrist

ReducFit = lm(BodyFat ~ catAge + Abdomen + BMI + Neck + 
                 Thigh + Hip + Forearm + Wrist, data = df2)
summary(ReducFit)
residualPlots(ReducFit3, type = "rstudent")
par(mfrow = c(2,1))
hist(ReducFit$resid, main = "Normality Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit_t_i = rstudent(ReducFit) 
qqnorm(ReducFit_t_i)                                         
qqline(ReducFit_t_i, col = 2, lwd = 2)
vif(ReducFit)
summary(influence.measures(ReducFit))
AIC(ReducFit, k = 2)
BIC(ReducFit)
ols_mallows_cp(ReducFit, fit0)
PRESS(ReducFit)
model_fit_stats(ReducFit)
summary(influence.measures(ReducFit))
residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot
anova(ReducFit)

#2nd Model Regressors are Abdomen, catAge, Neck, Thigh, Hip, Forearm, Wrist
ReducFit2 = lm(BodyFat ~ catAge + Abdomen + Neck + 
              Thigh + Hip + Forearm + Wrist, data = df2)
summary(ReducFit2)
residualPlots(ReducFit2, type = "rstudent")
par(mfrow = c(2,1))
hist(ReducFit2$resid, main = "Normalilty Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit2_t_i = rstudent(ReducFit2) 
qqnorm(ReducFit2_t_i)                                         
qqline(ReducFit2_t_i, col = 2, lwd = 2)
vif(ReducFit2)
summary(influence.measures(ReducFit2))
AIC(ReducFit2, k = 2)
BIC(ReducFit2)
ols_mallows_cp(ReducFit2, fit0)
PRESS(ReducFit2)
summary(influence.measures(ReducFit2))
residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit2) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot

#3rd Model Regressors are Abdomen, Wrist, catAge, BMI, Chest, Hip, Neck, Forearm and Thigh
ReducFit3 = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest + Hip + Neck + Forearm + Thigh,
               data = df3_ob39rem) 
summary(ReducFit3)
vif(ReducFit3)
summary(influence.measures(ReducFit3))
residualPlots(ReducFit3, type = "rstudent")
par(mfrow = c(2,1))
hist(ReducFit3$resid, main = "Normalilty Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit3_t_i = rstudent(ReducFit3) 
qqnorm(ReducFit3_t_i)                                         
qqline(ReducFit3_t_i, col = 2, lwd = 2)
residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit3) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot
AIC(ReducFit3, k = 2)
BIC(ReducFit3)
ols_mallows_cp(ReducFit3, alt_fit0)
PRESS(ReducFit3)


#4th Model Regressors are catAge, BMI, Chest, Abdomen, Wrist 
ReducFit4 = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest, data = df3_ob39rem)
summary(ReducFit4)
vif(ReducFit4)
par(mfrow = c(2,1))
hist(ReducFit4$resid, main = "Normalilty Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit4_t_i = rstudent(ReducFit4) 
qqnorm(ReducFit4_t_i)                                         
qqline(ReducFit4_t_i, col = 2, lwd = 2)
summary(influence.measures(ReducFit4))
residualPlots(ReducFit4, type = "rstudent")

residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit4) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot
AIC(ReducFit4, k = 2)
BIC(ReducFit4)
ols_mallows_cp(ReducFit4, alt_fit0)
PRESS(ReducFit4)

#5th Model Regressors are catAge, BMI, Chest, Abdomen, Hip and Wrist 
ReducFit5 = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest + Hip, data = df3_ob39rem)
summary(ReducFit5)
vif(ReducFit5)
residualPlots(ReducFit5, type = "rstudent")
par(mfrow = c(2,1))
hist(ReducFit5$resid, main = "Normalilty Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit5_t_i = rstudent(ReducFit5) 
qqnorm(ReducFit5_t_i)                                         
qqline(ReducFit5_t_i, col = 2, lwd = 2)
summary(influence.measures(ReducFit5))

residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit5) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot
AIC(ReducFit5, k = 2)
BIC(ReducFit5)
ols_mallows_cp(ReducFit5, alt_fit0)
PRESS(ReducFit5)

#6th Model Regressors are catAge, BMI, Neck, Chest, Abdomen, Hip, Forearm and Wrist 
ReducFit6 = lm(BodyFat ~ Abdomen + Wrist + catAge + BMI + Chest + Hip + Forearm + Neck, data = df3_ob39rem)
summary(ReducFit6)
vif(ReducFit6)
summary(influence.measures(ReducFit6))
residualPlots(ReducFit6, type = "rstudent")
par(mfrow = c(2,1))
hist(ReducFit6$resid, main = "Normalilty Plot for Body Fat Percentage", xlab ="Fitted Values")
ReducFit6_t_i = rstudent(ReducFit6) 
qqnorm(ReducFit6_t_i)                                         
qqline(ReducFit6_t_i, col = 2, lwd = 2)
residPlot = ggplot(aes(x = .fitted, y = .resid), data = ReducFit6) +
  geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
residPlot
AIC(ReducFit6, k = 2)
BIC(ReducFit6)
ols_mallows_cp(ReducFit6, alt_fit0)
PRESS(ReducFit6)

#UNUSED CODE #################################################################

# # Original Model
# #Response Variables: BodyFat
# #Explanatory Variables: Age, Weight, Height, Neck, Chest, 
# #Thigh, Abdomen, Hip, Knee, Ankle, Biceps, Forearm, Wrist
# 
# fit0 = lm(BodyFat ~ ., data = df1)          #fast way to get model
# summary(fit0)
# residualPlots(fit0, type = "rstudent")
# 
# par(mfrow = c(2,1))
# hist(fit0$resid, main = "Normalilty Plot for Body Fat Percentage")
# t_i = rstudent(fit0) 
# qqnorm(t_i)                                         # QQ plot
# qqline(t_i, col = 2, lwd = 2)
# skewness(fit0$residuals)
# #print(influence.measures(fitModel2))  #Prints all the measures for influence
# summary(influence.measures(fit0))
# cor(Data1)
corrplot(cor(df0), type = "upper",
         order = "original",
         method = "shade",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         diag = FALSE)
coefplot(fit0)
# head(fortify(fit0))
# residPlot = ggplot(aes(x = .fitted, y = .resid), data = fit0) 
# + geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
# vif(fit0)



#QQ plot looks pretty good
#Residual plot does show some pattern, a parabola type shape, may require transformation


# Transformed Model: square root transformation on BodyFat and Density
# sqrt_Density = sqrt(Data$Density)
# sqrt_BodyFat = sqrt(Data$BodyFat)
# 
# sqrt_fit0 = lm(BodyFat ~ sqrt_Density + Age + Weight + Height + Neck + Chest + 
#             Thigh + Abdomen + Hip + Knee + Ankle + Biceps
#             + Forearm + Wrist)
# summary(sqrt_fit0)
# residualPlots(sqrt_fit0, type = "rstudent")
# sqrt_t_i = rstudent(sqrt_fit0) 
# qqnorm(sqrt_t_i)   
# qqline(sqrt_t_i, col = 2, lwd = 2)
# #Square root has defintely made the entire model worse, although we can see that some regressors have become more significant.
# 
# # Transformed Model: reciprocal transformation on Density
# recip_Density = 1/(Data$Density)
# 
# recip_fit0 = lm(BodyFat ~ recip_Density + Age + Weight + Height + Neck + Chest + 
#                 Thigh + Abdomen + Hip + Knee + Ankle + Biceps
#                + Forearm + Wrist)
# 
# summary(recip_fit0)
# residualPlots(recip2_fit0, type = "rstudent")
# 
# 
# 
# # Transformed Model: reciprocal squared transformation on Density
# recip2_Density = 1/(Data$Density)^2
# 
# recip2_fit0 = lm(BodyFat ~ recip2_Density + Age + Weight + Height + Neck + Chest + 
#                 Thigh + Abdomen + Hip + Knee + Ankle + Biceps
#                 + Forearm + Wrist)
# 
# summary(recip2_fit0)
# residualPlots(recip2_fit0, type = "rstudent")
# 
# 
# # Transformed Model: reciprocal transformation on Density and Wrist to the power of 4
# recip4_Density = 1/(Data$Density)^4
# recip4_Wrist   = 1/(Data$Wrist)^4
# 
# recip4_fit0 = lm(BodyFat ~ recip4_Density + Age + Weight + Height + Neck + Chest + 
#                   Thigh + Abdomen + Hip + Knee + Ankle + Biceps
#                  + Forearm + recip4_Wrist)
# 
# summary(recip4_fit0)
# residualPlots(recip4_fit0, type = "rstudent")
# recip4_t_i = rstudent(recip3_fit0) 
# qqnorm(recip4_t_i)                                         # QQ plot
# qqline(t_i, col = 2, lwd = 2)











#Backward Selection based on AIC
# backward.fit1 =  step(fit0, direction='backward')
# summary(backward.fit1)
# backward.fit1$anova
# backward.fit1$coefficients
# vif(backward.fit1)
# residualPlots(backward.fit1, type = "rstudent")
# backward_t_i = rstudent(backward.fit1) 
# qqnorm(backward_t_i)                                         
# qqline(backward_t_i, col = 2, lwd = 2)
# 
# residPlot3 = ggplot(aes(x = .fitted, y = .resid), data = backward.fit1) +
#             geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")

#Backward Selection based on BIC
# backward.fit2 =  step(fit0, direction='backward', k = log(249))
# summary(backward.fit2)
# backward.fit2$anova
# backward.fit2$coefficients
# vif(backward.fit2)
# residualPlots(backward.fit2, type = "rstudent")
# backward_t_i = rstudent(backward.fit1) 
# qqnorm(backward_t_i)                                         
# qqline(backward_t_i, col = 2, lwd = 2)
# 
# residPlot3 = ggplot(aes(x = .fitted, y = .resid), data = backward.fit2) +
#   geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")



# 
# #Stepwise Regression based on AIC
# 
# step.fit1 = stepAIC(newFit0, direction = "both", scope = formula(newFit0))
# summary(step.fit1)
# step.fit1$anova
# vif(step.fit1)
# residualPlots(backward.fit1, type = "rstudent")
# step_t_i = rstudent(step.fit1) 
# qqnorm(step_t_i)                                         
# qqline(step_t_i, col = 2, lwd = 2)
# residPlot4 = ggplot(aes(x = .fitted, y = .resid), data = step.fit1) +
#   geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
# 
# #Stepwise Regression based on BIC
# fit2 = fit0
# step.fit2 = stepAIC(fit1, direction = "both", scope = formula(fit0), k = log(249))
# summary(step.fit2)
# step.fit2$anova
# vif(step.fit2)
# residualPlots(backward.fit2, type = "rstudent")
# step_t_i = rstudent(step.fit2) 
# qqnorm(step_t_i)                                         
# qqline(step_t_i, col = 2, lwd = 2)
# residPlot4 = ggplot(aes(x = .fitted, y = .resid), data = step.fit2) +
#   geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")





#Forward Selection
# intercept_only = lm(BodyFat ~ 1, data = df2)
# forward.fit1 = stepAIC(intercept_only, direction = 'forward', scope= formula(newFit0))
# summary(forward.fit1)
# forward.fit1$anova
# forward.fit1$coefficients
# vif(forward.fit1)
# residualPlots(forward.fit1, type = "rstudent")
# forward_t_i = rstudent(forward.fit1) 
# qqnorm(forward_t_i)                                         
# qqline(forward_t_i, col = 2, lwd = 2)
# residPlot2 = ggplot(aes(x = .fitted, y = .resid), data = forward.fit1) +
#   geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")
# 
# 
# #Forward Selection based on BIC
# forward.fit2 = stepAIC(intercept_only, direction = 'forward', scope= formula(fit0), k = log(249))
# summary(forward.fit2)
# forward.fit2$anova
# forward.fit2$coefficients
# vif(forward.fit2)
# residualPlots(forward.fit2, type = "rstudent")
# forward_t_i = rstudent(forward.fit2) 
# qqnorm(forward_t_i)                                         
# qqline(forward_t_i, col = 2, lwd = 2)
# residPlot2 = ggplot(aes(x = .fitted, y = .resid), data = forward.fit2) +
#   geom_point() + geom_hline(yintercept = 0) + labs(x = "Fitted Values", y = "Residuals")



# 
# #1st Model Regressors are Abdomen Weight Wrist Forearm Neck Age Thigh and Hip
# ReducFit = lm(BodyFat ~ Abdomen + Weight + Wrist + Forearm + Neck + Age + Thigh + Hip)
# summary(ReducFit)
# anova(ReducFit)
# residualPlots(ReducFit, type = "rstudent")
# par(mfrow = c(2,1))
# hist(ReducFit$resid, main = "Normalilty Plot for Body Fat Percentage")
# ReducFit_t_i = rstudent(bwdModFit) 
# qqnorm(ReducFit_t_i)                                         
# qqline(ReducFit_t_i, col = 2, lwd = 2)
# skewness(ReducFit$residuals)
# vif(ReducFit)
# AIC(ReducFit, k = 2)
# BIC(ReducFit)
# ols_mallows_cp(ReducFit, fit0)
# model_fit_stats(ReducFit)
# summary(influence.measures(ReducFit))
# 
# 
# #2nd Model Regressors are Weight Abdomen Forearm and Wrist
# ReducFit1 = lm(BodyFat ~ Weight + Abdomen + Forearm + Wrist)
# summary(ReducFit1)
# anova(ReducFit1)
# residualPlots(ReducFit1, type = "rstudent")
# par(mfrow = c(2,1))
# hist(ReducFit1$resid, main = "Normalilty Plot for Body Fat Percentage")
# ReducFit1_t_i = rstudent(bwdModFit) 
# qqnorm(ReducFit1_t_i)                                         
# qqline(ReducFit1_t_i, col = 2, lwd = 2)
# skewness(ReducFit1$residuals)
# vif(ReducFit1)
# AIC(ReducFit1, k = 2)
# BIC(ReducFit1)
# ols_mallows_cp(ReducFit1, fit0)
# PRESS(ReducFit)
# model_fit_stats(ReducFit1)
# summary(influence.measures(ReducFit1))


# #3rd Model Regressors are Age Weight Neck Thgigh Abdomen Forearm and Wrist
# ReducFit2 = lm(BodyFat ~ Age + Weight + Neck + Thigh + Abdomen + Forearm + Wrist)
# summary(ReducFit2)
# anova(ReducFit2)
# residualPlots(ReducFit2, type = "rstudent")
# par(mfrow = c(2,1))
# hist(ReducFit2$resid, main = "Normalilty Plot for Body Fat Percentage")
# ReducFit2_t_i = rstudent(bwdModFit) 
# qqnorm(ReducFit2_t_i)                                         
# qqline(ReducFit2_t_i, col = 2, lwd = 2)
# skewness(ReducFit2$residuals)
# vif(ReducFit2)
# AIC(ReducFit2, k = 2)
# BIC(ReducFit2)
# ols_mallows_cp(ReducFit2, fit0)
# model_fit_stats(ReducFit2)
# summary(influence.measures(ReducFit2))

# #Interaction
# interactFit0 = lm(BodyFat ~ catAge + BMI + Neck + Chest + Thigh + Abdomen + 
#                     Hip + Knee + Ankle + Biceps + Forearm + Wrist + BMI*catAge, data = df2)
# 
# summary(interactFit0)
# interaction.plot(BMI, catAge, Data1$BodyFat,
#                  xlab = "Body Mass Index",
#                  ylab = "Body Fat (%)",
#                  main = "Interaction between Different Age Cohorts and Body Fat (%)",
#                  trace.label = "Age Group", type = "b", col=c("red", "green", "black"), 
#                  pch= c(19,17), fixed = TRUE)