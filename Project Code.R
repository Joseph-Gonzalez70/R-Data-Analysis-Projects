#Data set 1 bone marrow
#Read in the data
data1= read.table("T7-1.DAT", header=FALSE)
colnames(data1)=c("Total_Dwelling_Size", "Assessed_Value", "Selling_Price")
attach(data1)

#Summary of data
sapply(data1,class)
summary(data1)
sapply(data1,var)
par(mfrow=c(1,3))
xlabels=c("Hundred Square Feet", "Thouand Dollars", "Thousand Dollars")
#Histogram
for (i in 1:3){
  hist(data1[,i], main=paste("Histogram of", names(data1[i])),xlab=xlabels[i])
}
par(mfrow=c(2,2))
#Boxplot
for (i in 1:3){
  boxplot(data1[,i], main=paste("Boxplot of", names(data1[i])),ylab=xlabels[i])
}
#Histogram for residuals
fit1=lm(Selling_Price~ Total_Dwelling_Size+Assessed_Value, data = data1)
fit2=lm(Selling_Price~ Total_Dwelling_Size, data = data1)
anova(fit1)
par(mfrow=c(1,3))
hist(fit1$residuals, main="Histogram of Residuals")
plot(fit1, which=1)
plot(fit1, which=2)

-----------------------------------------------------------------------------------------------
#Tests
n <- length(Selling_Price)
Z <- cbind(rep(1,n),as.matrix(data1[,1:2]))
Z
sapply(Z,class)
r <- dim(Z)[2]-1

#Least square estimates
beta_hat <- solve(t(Z)%*%Z)%*%t(Z)%*%Selling_Price
beta_hat

#Sums of squares
#SSTO
Yadjust=data1$Selling_Price-mean(data1$Selling_Price)*rep(1,n)
SSTO=sum(Yadjust^2)

#ESS
Yadjusted2=fit1$fitted.values-mean(data1$Selling_Price)*rep(1,n)
SSRtotal=sum(Yadjusted2^2)
SSRtotal

#RSS
Rres1=0
sum(fit1$residuals^2)

#R squared value
R_square <- 1 - sum((Selling_Price - Z%*%beta_hat)^2)/sum((Selling_Price-mean(Selling_Price))^2)
R_square

#Sigma_hat_square
sigma_hat_square <- sum((Selling_Price - Z%*%beta_hat)^2)/(n-r-1)
sigma_hat_square

#Estimated covariance of hat{beta}
cov_B = sigma_hat_square * solve(t(Z)%*%Z)
cov_B

#Confidence interval for beta_1
alpha=0.05
j <- 2
cat('[',
    beta_hat[j+1] - qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/2, n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')


#Confidence region based simultaneous confidence intervals 

j <- 0
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

j <- 1
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

j <- 2
cat('[',
    beta_hat[j+1] - sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + sqrt((r+1)*qf(1-alpha, r+1, n-r-1))*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

#Bonferroni correction based simultaneous confidence intervals

j <- 0
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

j <- 1
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')


j <- 2
cat('[',
    beta_hat[j+1] - qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ',',
    beta_hat[j+1] + qt(1-alpha/(2*(r+1)), n-r-1)*sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1]),
    ']')

#F-test
#H_0: beta_1 = beta_2 = 0

C <- matrix(c(0,0,0,0,1,1),2,3)
C

df_1 <- qr(C)$rank # df_1: rank of matrix R

f_stat <- (t(C%*%beta_hat)%*%solve(C%*%solve(t(Z)%*%Z)%*%t(C))%*%(C%*%beta_hat)/df_1)/sigma_hat_square
f_stat

cval_f <- qf(1-alpha, 2, n-r-1)
cval_f

anova(fit1,fit2)
qf(1-0.05,1,17)
#t-test for single coefficient
#H_0: beta_j = 0, H_a: beta_j != 0

j <- 1
t_stat <- (beta_hat[j+1] - 0)/sqrt(sigma_hat_square * solve(t(Z)%*%Z)[j+1,j+1])
t_stat

alpha <- 0.05
cval_t <- qt(1-alpha/2, n-r-1)
cval_t

#Confidence interval for z_0^T beta

z_0 <- c(1, mean(Z[,2]), mean(Z[,3]))

cat('[',
    z_0%*%beta_hat - sqrt(sigma_hat_square)*sqrt(t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ',',
    z_0%*%beta_hat + sqrt(sigma_hat_square)*sqrt(t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ']')

#Prediction interval for Y_0 = z_0^T beta + epsilon_0

cat('[',
    z_0%*%beta_hat - sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ',',
    z_0%*%beta_hat + sqrt(sigma_hat_square)*sqrt(1+t(z_0)%*%solve(t(Z)%*%Z)%*%z_0)*qt(1-alpha/2, n-r-1),
    ']')

#Confidence Region for (beta_1, beta_2)^T

center <- beta_hat[2:3]
es<-eigen(C%*%solve(t(Z)%*%Z)%*%t(C))
e1<-es$vec %*% diag(sqrt(es$val))
r1<-sqrt(df_1*cval_f*sigma_hat_square)
theta<-seq(0,2*pi,len=250)
v1<-cbind(r1*cos(theta), r1*sin(theta))
pts<-t(center - (e1%*%t(v1)))
plot(pts,type="l",main="Confidence Region for (beta_1, beta_2)^T",xlab="beta_1",ylab="beta_2",asp=1,
     xlim = c(-1,5),ylim=c(-1,1))
segments(0,center[2],center[1],center[2],lty=2) # highlight the center
segments(center[1],0,center[1],center[2],lty=2)
arrows(-1,0,5,0)
arrows(0,-2,0,2)

th2<-c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
v2<-cbind(r1*cos(th2), r1*sin(th2))
pts2<-t(center-(e1%*%t(v2)))
segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)

##########################################################################################
------------------------------------------------------------------------------------------------------
#Dataset 2
#Data
data2= read.table("T6-9.DAT", header=FALSE)
colnames(data2)=c("Length", "Width", "Height", "Gender")
attach(data2)

#Boxplot
for (i in 1:3){
  boxplot(data2[,i]~data2$Gender, main=paste("Boxplot of", names(data2[i])),ylab=xlabels[i])
}

#############Analysis

#### two-sample Hotelling's T2 test  ---------------------

males=which(data2$Gender=="male")
females=which(data2$Gender=="female")
Maledata <- data2[males,1:3]
Femaledata <- data2[females,1:3]
sapply(data2,class)
summary(Maledata)
summary(Femaledata)
sapply(Maledata[,1:3],var)
sapply(Maledata[,1:3],sd)
sapply(Femaledata[,1:3],var)
sapply(Femaledata[,1:3],sd)
par(mfrow=c(1,3))
xlabels=c("Length", "Width", "Height")

#Now, we perform the two-sample Hotelling T^2-test
n<-c(24,24)
p<-3
xmean1<-colMeans(Maledata)
xmean2<-colMeans(Femaledata)
d<-xmean1-xmean2
S1<-var(Maledata)
S2<-var(Femaledata)
Sp<-((n[1]-1)*S1+(n[2]-1)*S2)/(sum(n)-2)
t2 <- t(d)%*%solve(sum(1/n)*Sp)%*%d
t2

alpha<-0.05
cval <- (sum(n)-2)*p/(sum(n)-p-1)*qf(1-alpha,p,sum(n)-p-1)
cval
sum(n)-p-1

#To check the significant components
#Simultaneous confidence intervals
wd<-sqrt(cval*diag(Sp)*sum(1/n))
Cis<-cbind(d-wd,d+wd)

cat("95% simultaneous confidence interval","\n")
Cis

#Bonferroni simultaneous confidence intervals
wd.b<- qt(1-alpha/(2*p),n[1]+n[2]-2) *sqrt(diag(Sp)*sum(1/n))
Cis.b<-cbind(d-wd.b,d+wd.b)
cat("95% Bonferroni simultaneous confidence interval","\n")
Cis.b

-----------------------------------------------------------------------------------
#LDA Analysis on Dataset 2
library(rrcov)

par(mfrow=c(1,3), mar=c(4,4,2,1))
plot(data2$Length,data2$Width,xlab="Length",ylab="Width",
     pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Width Vs. Length")
legend("topright",legend=c("Females","Males"),pch=c(18,20),col=c(2,4),cex=1)

plot(data2$Length,data2$Height,xlab="Length",ylab="Height",
     pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Height Vs. Length")
legend("topright",legend=c("Females","Males"),pch=c(18,20),col=c(2,4),cex=1)

plot(data2$Width,data2$Height,xlab="Width",ylab="Height",
     pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Height Vs.Width")
legend("topright",legend=c("Females","Males"),pch=c(18,20),col=c(2,4),cex=1)


##Method 2: use function LDA in MASS package:

#Width Vs Length----------------------------------------------------------
library(MASS)
lda.obj <- lda(Gender~Length+ Width,data=data2,prior=c(1,1)/2)
plda <- predict(object=lda.obj,newdata=data2)

#Confusion matrix
table(data2[,4],plda$class)

#Expected actual error rate
n <- dim(data2)[1]
n_M <- 0
for (i in 1:n){
  lda.obj <- lda(Gender~ Length + Width,data=data2[-c(i),],prior=c(1,1)/2)
  plda <- predict(object=lda.obj,data2[c(i),])
  n_M <- n_M + (plda$class != data2[c(i),]$Gender)
}
n_M/n

#Plot the decision line
gmean <- lda.obj$prior %*% lda.obj$means
const <- as.numeric(gmean %*%lda.obj$scaling)
slope <- - lda.obj$scaling[1] / lda.obj$scaling[2]
intercept <- const / lda.obj$scaling[2]

#Plot decision boundary
plot(data2[,1:2],pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Width Vs. Length")
abline(intercept, slope)
legend("topright",legend=c("Female","Male"),pch=c(18,20),col=c(2,4))

#Height V Length----------------------------------------------------------
lda.obj <- lda(Gender~Length+ Height,data=data2,prior=c(1,1)/2)
plda <- predict(object=lda.obj,newdata=data2)

?lda
#Confusion matrix
table(data2[,4],plda$class)

#Expected Actual Error Rate
n <- dim(data2)[1]
n_M <- 0
for (i in 1:n){
  lda.obj <- lda(Gender~ Length + Height,data=data2[-c(i),],prior=c(1,1)/2)
  plda <- predict(object=lda.obj,data2[c(i),])
  n_M <- n_M + (plda$class != data2[c(i),]$Gender)
}
n_M/n


#Plot the decision line
gmean <- lda.obj$prior %*% lda.obj$means
const <- as.numeric(gmean %*%lda.obj$scaling)
slope <- - lda.obj$scaling[1] / lda.obj$scaling[2]
intercept <- const / lda.obj$scaling[2]

#Plot decision boundary
plot(data2[,c(1,3)],pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Height Vs. Length")
abline(intercept, slope)
legend("topright",legend=c("Female","Male"),pch=c(18,20),col=c(2,4))


#Height V. Width-------------------------------------------------------
?lda
lda.obj <- lda(Gender~Width+ Height,data=data2,prior=c(1,1)/2)
plda <- predict(object=lda.obj,newdata=data2)

# Confusion matrix
table(data2[,4],plda$class)

#Expected Actual Error Rate
n <- dim(data2)[1]
n_M <- 0
for (i in 1:n){
  lda.obj <- lda(Gender~ Width + Height,data=data2[-c(i),],prior=c(1,1)/2)
  plda <- predict(object=lda.obj,data2[c(i),])
  n_M <- n_M + (plda$class != data2[c(i),]$Gender)
}
n_M/n

#Plot the decision line
gmean <- lda.obj$prior %*% lda.obj$means
const <- as.numeric(gmean %*%lda.obj$scaling)
slope <- - lda.obj$scaling[1] / lda.obj$scaling[2]
intercept <- const / lda.obj$scaling[2]

#Plot decision boundary
plot(data2[,2:3],pch=rep(c(18,20),each=24),col=rep(c(2,4),each=24),main="Height Vs.Width")
abline(intercept, slope)
legend("topright",legend=c("Female","Male"),pch=c(18,20),col=c(2,4))


--------------------------------------------------------------------------------------------------------
  ##########################################################################################
------------------------------------------------------------------------------------------------------
#dataset 3 PCA
data3= read.table("T8-5.DAT", header=FALSE)
colnames(data3)=c("Tot_population", "Profess_deg", "employed_16", "Gov_employ", "Med_home_value")
attach(data3)
sapply(data3,class)
summary(data3)
sapply(data3,var)
sapply(data3,sd)

#Correlation matrix
census.pc <- princomp(data3, cor=TRUE)

summary(census.pc, loadings = TRUE)

#Showing the eigenvalues of the correlation matrix:
(census.pc$sdev)^2

#A scree plot:
plot(1:(length(census.pc$sdev)),  (census.pc$sdev)^2, type='b', 
     main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")

#Plotting the PC scores for the sample data in the space of the first two principal components:
par(pty="s")
plot(census.pc$scores[,1], census.pc$scores[,2], 
     xlab="PC 1", ylab="PC 2", type ='n', lwd=2, main="Principle Component Scores")





