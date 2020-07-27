## partial code credit to Weiping Zhang(USTC)
#### two-sample Hotelling's T2 test  -------

Setosa.sepal <- iris[iris$Species == "setosa",-c(3,4,5)]
Versicolor.sepal <- iris[iris$Species == "versicolor",-c(3,4,5)]

# now we perform the two-sample Hotelling T^2-test
n<-c(50,50)
p<-2
xmean1<-colMeans(Setosa.sepal)
xmean2<-colMeans(Versicolor.sepal)
d<-xmean1-xmean2
S1<-var(Setosa.sepal)
S2<-var(Versicolor.sepal)
Sp<-((n[1]-1)*S1+(n[2]-1)*S2)/(sum(n)-2)
t2 <- t(d)%*%solve(sum(1/n)*Sp)%*%d
t2

alpha<-0.05
cval <- (sum(n)-2)*p/(sum(n)-p-1)*qf(1-alpha,p,sum(n)-p-1)
cval



# Confidence Region
es<-eigen(sum(1/n)*Sp)
e1<-es$vec %*% diag(sqrt(es$val))
r1<-sqrt(cval)
theta<-seq(0,2*pi,len=250)
v1<-cbind(r1*cos(theta), r1*sin(theta))
pts<-t(d-(e1%*%t(v1)))
plot(pts,type="l",main="Confidence Region for Bivariate Normal",xlab="Sepal.length",ylab="Sepal.width",asp=1)
segments(0,d[2],d[1],d[2],lty=2) # highlight the center
segments(d[1],0,d[1],d[2],lty=2)

th2<-c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
v2<-cbind(r1*cos(th2), r1*sin(th2))
pts2<-t(d-(e1%*%t(v2)))
segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)



# since we reject the null, we use the simultaneous confidence intervals
# to check the significant components

# simultaneous confidence intervals
wd<-sqrt(cval*diag(Sp)*sum(1/n))
Cis<-cbind(d-wd,d+wd)

cat("95% simultaneous confidence interval","\n")
Cis

#Bonferroni simultaneous confidence intervals
wd.b<- qt(1-alpha/(2*p),n[1]+n[2]-2) *sqrt(diag(Sp)*sum(1/n))
Cis.b<-cbind(d-wd.b,d+wd.b)
cat("95% Bonferroni simultaneous confidence interval","\n")
Cis.b

# both component-wise simultaneous confidence intervals do not contain 0, so they have significant differences. 