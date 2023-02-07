{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 CourierNewPSMT;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs29\fsmilli14667 \cf2 \expnd0\expndtw0\kerning0
rho=0.3\
d = 10^6\
z<-rnorm(d)\
u1<-rnorm(d)\
u2<-rnorm(d)\
u3<-rnorm(d)\
x1=sqrt(rho)*z+sqrt(1-rho)*u1\
x2=sqrt(rho)*z+sqrt(1-rho)*u2\
x3=sqrt(rho)*z+sqrt(1-rho)*u3\
data <- data.frame(x1,x2,x3)\
\
mean(x1);mean(x2);mean(x3)\
\pard\pardeftab720\sa240\partightenfactor0
\cf2 cov(data) \
\pard\pardeftab720\partightenfactor0
\cf2 plot(x1,x2)\
plot(x1,x3)\
plot(x2,x3)\
}