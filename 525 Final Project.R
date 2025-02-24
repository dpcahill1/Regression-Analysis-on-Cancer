#packages
library(car)
library(ggplot2)
library(plotly)
library(interactions)
library(olsrr)

#data
data <- read.csv("~/Desktop/STAT 525/Data Sets/cancer.csv", header = T)
data$Grade <-ifelse(data$ClassGrade == 'UDH', "UDH", 'nUDH')

#misc values
n <- nrow(data)
p <- 10

upper <- 3*p/n
lower <- 2*p/n

#initial model and assumption
model1 <- lm(Mean_Area ~ Mean_mean_R + Mean_mean_G + Mean_mean_B + Mean_mean_HE + Mean_mean_BR + Mean_mean_Lab + Mean_mean_HSV * Grade, data = data)
outdata1 <- fortify(model1)

pairs(~ Mean_Area + Mean_mean_R + Mean_mean_G + Mean_mean_B + Mean_mean_HSV + 
        Mean_mean_HE + Mean_mean_BR + Mean_mean_Lab, data=data, main='Scatterplot Matrix')

ggplot(outdata1, aes(x=.fitted, y=.stdresid)) +
  geom_point(shape=16, size=3) +
  labs(x = "Predicted Mean Area",y="Studentized Residuals",
       title="Residual against Predicted Response")+
  theme_bw()+
  geom_hline(yintercept=0)+
  theme(plot.title =
          element_text(hjust=0.5, size=rel(1.6)))+
  theme(axis.title.y = element_text(size = rel(1.4)))+
  theme(axis.title.x = element_text(size = rel(1.4)))+
  theme(axis.text.x = element_text(size = rel(1.6)))+
  theme(axis.text.y = element_text(size = rel(1.6)))

outdata1$yHatCategory <- ifelse(
  outdata1$.fitted < median(outdata1$.fitted), c("group1"), c("group2"))
leveneTest(.resid ~ yHatCategory, data=outdata1)

par(mfrow=c(1,2))
hist(outdata1$.resid,
     main = 'Histogram of Residuals', xlab='Residuals')
qqnorm(outdata1$.resid, main = "Normal Q-Q Plot",
       xlab = "Theoretical Normal Quantiles",
       ylab = "Sample Residuals", pch = 16)
qqline(outdata1$.resid)

shapiro.test(outdata1$.resid)

levplot1 <- ggplot(outdata1,aes(x=.hat,y=.stdresid,size=.cooksd)) +
  geom_point(shape=16) +
  labs(x='Leverage Values', y = 'Studentized Residuals',
       title = 'Leverage and Influence Plot')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0.1198,color='red')+
  annotate("text", x = 0.25, y = -2,
           label='leverage cutoff = 0.1198')+
  geom_vline(xintercept = 0.1796,color='red')+
  annotate("text", x = 0.4, y = -2,
           label='leverage cutoff = 0.1796')+
  theme(plot.title = element_text(hjust=0.5,size=rel(1.6)))+
  theme(axis.title.y = element_text(size = rel(1.4)))+
  theme(axis.title.x = element_text(size = rel(1.4)))+
  theme(axis.text.x = element_text(size = rel(1.6)))+
  theme(axis.text.y = element_text(size = rel(1.6)))
interactiveplot1 <-ggplotly(levplot1)
interactiveplot1

intmodel2 <- lm(Mean_Area ~ Mean_mean_HSV*Grade, data = data)
summary(intmodel2)
interact_plot(intmodel2, pred = 'Mean_mean_HSV', modx = 'Grade',
              plot.points=TRUE, interval=FALSE, x.label='Average Lightness Across Nuclei on HSV Scale',
              y.label='Mean Area', main.title='Mean Area vs. Average Lightness Across Nuclei (by Class Grade)',
              legend.main='', colors = 'Blues', line.thickness=1)

#long transformed model and assumptions
data$MeanAreaTransformed <- log(data$Mean_Area)
model <- lm(MeanAreaTransformed ~ Mean_mean_R + Mean_mean_G + Mean_mean_B + Mean_mean_HE + Mean_mean_BR + Mean_mean_Lab + Mean_mean_HSV * Grade, data = data)
outdata <- fortify(model)

pairs(~ MeanAreaTransformed + Mean_mean_R + Mean_mean_G + Mean_mean_B + Mean_mean_HSV + 
        Mean_mean_HE + Mean_mean_BR + Mean_mean_Lab, data=data, main='Scatterplot Matrix')

ggplot(outdata, aes(x=.fitted, y=.stdresid)) +
  geom_point(shape=16, size=3) +
  labs(x = "Predicted Mean Area",y="Studentized Residuals",
       title="Residual against Predicted Response")+
  theme_bw()+
  geom_hline(yintercept=0)+
  theme(plot.title =
          element_text(hjust=0.5, size=rel(1.6)))+
  theme(axis.title.y = element_text(size = rel(1.4)))+
  theme(axis.title.x = element_text(size = rel(1.4)))+
  theme(axis.text.x = element_text(size = rel(1.6)))+
  theme(axis.text.y = element_text(size = rel(1.6)))

outdata$yHatCategory <- ifelse(
  outdata$.fitted < median(outdata$.fitted), c("group1"), c("group2"))
leveneTest(.resid ~ yHatCategory, data=outdata)

par(mfrow=c(1,2))
hist(outdata$.resid,
     main = 'Histogram of Residuals', xlab='Residuals')
qqnorm(outdata$.resid, main = "Normal Q-Q Plot",
       xlab = "Theoretical Normal Quantiles",
       ylab = "Sample Residuals", pch = 16)
qqline(outdata$.resid)

shapiro.test(outdata$.resid)

levplot <- ggplot(outdata,aes(x=.hat,y=.stdresid,size=.cooksd)) +
  geom_point(shape=16) +
  labs(x='Leverage Values', y = 'Studentized Residuals',
       title = 'Leverage and Influence Plot')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0.1198,color='red')+
  annotate("text", x = 0.25, y = -2,
           label='leverage cutoff = 0.1198')+
  geom_vline(xintercept = 0.1796,color='red')+
  annotate("text", x = 0.4, y = -2,
           label='leverage cutoff = 0.1796')+
  theme(plot.title = element_text(hjust=0.5,size=rel(1.6)))+
  theme(axis.title.y = element_text(size = rel(1.4)))+
  theme(axis.title.x = element_text(size = rel(1.4)))+
  theme(axis.text.x = element_text(size = rel(1.6)))+
  theme(axis.text.y = element_text(size = rel(1.6)))
interactiveplot <-ggplotly(levplot)
interactiveplot

model2 <- lm(MeanAreaTransformed ~ Mean_mean_R + Mean_mean_G + Mean_mean_B +  Mean_mean_HE + Mean_mean_BR + Mean_mean_Lab , data = data)
ols_step_both_p(model2, p_val = 0.05, details = T)
seq_model<-lm(MeanAreaTransformed ~ Mean_mean_B +Mean_mean_BR+ Mean_mean_HSV + Grade + Mean_mean_HSV*Grade, data = data)
summary(seq_model)

#centered model
data$cMean_mean_HSV <- data$Mean_mean_HSV - mean(data$Mean_mean_HSV)
centermod <-lm(MeanAreaTransformed ~ Mean_mean_B + Mean_mean_BR + cMean_mean_HSV*Grade, data = data)
outdatacent <- fortify(centermod)

outlier <- which(outdatacent$.stdresid > 3 | outdatacent$.hat > 0.1617)
outdatacent <- outdatacent[-outlier,]

summary(centermod)

pairs(~ MeanAreaTransformed + Mean_mean_B + cMean_mean_HSV + Mean_mean_BR , data=data, main='Scatterplot Matrix')

ggplot(outdatacent, aes(x=.fitted, y=.stdresid)) +
  geom_point(shape=16, size=3) +
  labs(x = "Predicted Mean Area",y="Studentized Residuals",
       title="Residual against Predicted Response")+
  theme_bw()+
  geom_hline(yintercept=0)+
  theme(plot.title =
          element_text(hjust=0.5, size=rel(1.6)))+
  theme(axis.title.y = element_text(size = rel(1.4)))+
  theme(axis.title.x = element_text(size = rel(1.4)))+
  theme(axis.text.x = element_text(size = rel(1.6)))+
  theme(axis.text.y = element_text(size = rel(1.6)))

outdatacent$yHatCategory <- ifelse(
  outdatacent$.fitted < median(outdatacent$.fitted), c("group1"), c("group2"))
leveneTest(.resid ~ yHatCategory, data=outdatacent)

par(mfrow=c(1,2))
hist(outdatacent$.resid,
     main = 'Histogram of Residuals', xlab='Residuals')
qqnorm(outdatacent$.resid, main = "Normal Q-Q Plot",
       xlab = "Theoretical Normal Quantiles",
       ylab = "Sample Residuals", pch = 16)
qqline(outdatacent$.resid)

shapiro.test(outdatacent$.resid)

vif(centermod, type = "predictor")

centermod2 <-lm(MeanAreaTransformed ~ Mean_mean_BR + cMean_mean_HSV*Grade, data = data)
vif(centermod2, type = "predictor")

model5 <- lm(MeanAreaTransformed ~ Mean_mean_B + Mean_mean_BR + cMean_mean_HSV * Grade, data = data)
summary(model5)

