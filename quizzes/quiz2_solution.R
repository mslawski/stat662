rho=0.3
d = 10^6
z<-rnorm(d)
u1<-rnorm(d)
u2<-rnorm(d)
u3<-rnorm(d)
x1=sqrt(rho)*z+sqrt(1-rho)*u1
x2=sqrt(rho)*z+sqrt(1-rho)*u2
x3=sqrt(rho)*z+sqrt(1-rho)*u3
data <- data.frame(x1,x2,x3)

mean(x1);mean(x2);mean(x3)
cov(data) 
plot(x1,x2)
plot(x1,x3)
plot(x2,x3)
