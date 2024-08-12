ggplot(aes(x=d,color=b,fill=b),
       data=rbind(data.frame(b="large",d=rbeta(10000,1/10,1/10)+1),
                  data.frame(b="mediu",d=rbeta(10000,1/3*1/(1+3),1/3)+1),
                  data.frame(b="small",d=rbeta(10000,1/10*1/(1+1),1/10)+1)))+
  stat_density(alpha=0.2,position = position_dodge())

ggplot()+
  stat_function(fun = function(x) dbeta(x,1/10,1/10),n = 1000)+
  stat_function(fun = function(x) dbeta(x,1/10,1/3),color="red",n = 1000)+
  xlim(0.01,0.99)

ggplot(data.frame(x=1:32,y=0))+geom_line(aes(x=x,y=y))+
  stat_function(fun=function(x) dgeom(x,p=1/10),n=32)+
  stat_function(fun=function(x) dgeom(x,p=1/10),n=32)+
  stat_function(fun=function(x) dgeom(x,p=8/10),n=32,color="red")

x <- 21
y <- 5
z <- 2
ggplot(aes(x=d,color=b,fill=b),
       data=rbind(data.frame(b="many",d=rbeta(10000,1/x,1/x)),
                  data.frame(b="manyx",d=rbeta(10000,1/x,1/x)),
                  data.frame(b="few",d=rbeta(10000,1/z,1/z)),
                  data.frame(b="mid",d=rbeta(10000,1/y,1/y))))+
  stat_density(alpha=0.2,position = position_dodge())

ggplot(aes(x=d,color=b,fill=b),
       data=rbind(
         data.frame(b="rs=0.1, ts=0.1, l=1",d=sapply(rpois(10000,1), function (c)  (c+1)*(log2(rgeom(1, 0.1)+1) + (rgeom(1, 0.01)+1)))),
         data.frame(b="rs=0.3, ts=0.4, l=1",d=sapply(rpois(10000,1), function (c)  2*(c+1)*(log2(rgeom(1, 0.3)+1) + (rgeom(1, 0.04)+1)))),
         data.frame(b="rs=0.5, ts=0.7, l=2",d=sapply(rpois(10000,2), function (c)  (c+1)*(log2(rgeom(1, 0.5)+1) + (rgeom(1, 0.1)+1)))),
         data.frame(b="rs=0.6, ts=0.6, l=2",d=sapply(rpois(10000,2), function (c)  (c+1)*(log2(rgeom(1, 0.6)+1) + (rgeom(1, 0.2)+1)))),
         data.frame(b="rs=0.9, ts=0.85, l=3",d=sapply(rpois(10000,3), function (c)  (c+1)*(log2(rgeom(1, 0.2)+1) + (rgeom(1, 0.4)+1)))))
       )+stat_density(alpha=0.2, position = position_dodge()) + scale_x_continuous(trans = "log2")
       
       b <- 16
       y <- 21
       ggplot(aes(x=d,color=b,fill=b),
              data=rbind(
                data.frame(b="row=4, size=0.2, c=2",d=sapply(rbeta(10000,1/y,1/y), function (p)  2*(rbinom(1, ceiling(1/1*0.1*b), p)+1))),
                data.frame(b="row=2, size=0.2, c=5",d=sapply(rbeta(10000,1/y,1/y), function (p)  5*(rbinom(1, ceiling(1/2*0.2*b), p)+1))),
                data.frame(b="row=1, size=0.2, c=1",d=sapply(rbeta(10000,1/y,1/y), function (p)  1*(rbinom(1, ceiling(1/1*0.2*b), p)+1))),
                data.frame(b="row=8, size=0.5, c=3",d=sapply(rbeta(10000,1/y,1/y), function (p)  3*(rbinom(1, ceiling(1/8*0.5*b), p)+1))),
                data.frame(b="row=10, size=1.0, c=5",d=sapply(rbeta(10000,1/y,1/y), function (p)  4*(rbinom(1, ceiling(1/10*1.0*b), p)+1)))),
       )+
         stat_density(alpha=0.2, position = position_dodge())
       ggplot(aes(x=d,color=b,fill=b),
              data=rbind(
                data.frame(b="size=max*0.2, c=2",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 2*ceiling(0.2*b), p)+1)),
                data.frame(b="size=max*0.2, c=3",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 3*ceiling(0.2*b), p)+1)),
                data.frame(b="size=max*0.2, c=5",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 5*ceiling(0.2*b), p)+1)),
                data.frame(b="size=max*0.2, c=1",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 1*ceiling(0.2*b), p)+1)),
                data.frame(b="size=max*0.5, c=1",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 1*ceiling(0.5*b), p)+1)),
                data.frame(b="size=max*0.5, c=3",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 3*ceiling(0.5*b), p)+1)),
                data.frame(b="size=max*1.0, c=1",d=sapply(rbeta(10000,1/y,1/y), function (p) rbinom(1, 1*ceiling(1.0*b), p)+1)))
       )+ stat_density(alpha=0.2, position = position_dodge())
       
       ggplot(aes(x=d,color=b,fill=b),
              data=rbind(
                data.frame(b="size=max*0.2, c=2",d=sapply(rbeta(10000,1/2,1/2), function (p) 2*rnbinom(1, ceiling(0.2*b)+1, p))),
                data.frame(b="size=max*0.2, c=5",d=sapply(rbeta(10000,1/2,1/2), function (p) 5*rnbinom(1, ceiling(0.2*b)+1, p))),
                data.frame(b="size=max*0.2, c=1",d=sapply(rbeta(10000,1/2,1/2), function (p) 1*rnbinom(1, ceiling(0.2*b)+1, p))),
                data.frame(b="size=max*0.5, c=1",d=sapply(rbeta(10000,1/2,1/2), function (p) 1*rnbinom(1, ceiling(0.5*b)+1, p))),
                data.frame(b="size=max*0.5, c=3",d=sapply(rbeta(10000,1/2,1/2), function (p) 3*rnbinom(1, ceiling(0.5*b)+1, p))),
                data.frame(b="size=max*1.0, c=1",d=sapply(rbeta(10000,1/2,1/2), function (p) 1*rnbinom(1, ceiling(1.0*b)+1, p))))
       )+ stat_ecdf(alpha=1) + scale_x_continuous(trans = "log")
       
       ggplot(aes(x=d,color=b,fill=b),
              data=rbind(data.frame(b="large",d=rbeta(10000,1/y,1/y)),
                         data.frame(b="small",d=rbeta(10000,1/y,1/y))))+
         stat_density(alpha=0.2,position = position_dodge())
       