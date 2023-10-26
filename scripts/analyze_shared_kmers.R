require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

taxon_dist <- read_csv("../data/dist_wrt_matchrank.csv")

taxon_dist$rank <- factor(
  taxon_dist$rank,
  levels = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")
)

ggplot(taxon_dist %>% filter(rank !="kingdom" & rank != "root")) + aes(x=rank, y=d, color=rank) +
  geom_boxplot() +
  geom_point()
  

L=10^6
n=30

sharedkmers(30,0.01,n+1,10^6)/((L*n+1)-sharedkmers(30,0.01,n+1,10^6))

atleastone = function(k,d) 1-(1-d)^k
sharedkmers = function(k,d,n=n,L=L)  L*n - L*n * (atleastone(k,d)+((1-atleastone(k,d))*atleastone(k,d)^(n-1)))
sharedkmerprop = function(k,d,n=n,L=L) sharedkmers(k,d,n=n,L=L)/(n*L)

ggplot(aes(x=x),data=data.frame(x=1:40))+
  stat_function(aes(color=" 2% (species)"),fun=function(x) sharedkmerprop(30,0.01,x+1,10^6))+
  stat_function(aes(color=" 5% (species/genus)"),fun=function(x) sharedkmerprop(30,0.025,x+1,10^6))+
  stat_function(aes(color="10% (genus)"),fun=function(x) sharedkmerprop(30,0.05,x+1,10^6))+
  stat_function(aes(color="15% (genus/family)"),fun=function(x) sharedkmerprop(30,0.075,x+1,10^6))+
  stat_function(aes(color="20% (family)"),fun=function(x) sharedkmerprop(30,0.1,x+1,10^6))+
  stat_function(aes(color="25% (family/order)"),fun=function(x) sharedkmerprop(30,0.125,x+1,10^6))+
  stat_function(aes(color="30% (order)"),fun=function(x) sharedkmerprop(30,0.15,x+1,10^6))+
  stat_function(aes(color="35% (order/class)"),fun=function(x) sharedkmerprop(30,0.175,x+1,10^6))+
  stat_function(aes(color="40% (class)"),fun=function(x) sharedkmerprop(30,0.2,x+1,10^6))+
  theme_cowplot() +
  scale_color_brewer(palette = "Paired", direction = -1) +
  scale_x_continuous(trans = "log2")

L=10^6
n=10

sharedkmerprop(30, 0.01, 32, L)

atleastone = function(k,d) 1-(1-d)^k
sharedkmers = function(k,d,n=n,L=L) L*n*(atleastone(k,d)^(n)+((1-atleastone(k,d))*atleastone(k,d)^(n-1)))
sharedkmerprop = function(k,d,n=n,L=L) 1-sharedkmers(k,d,n=n,L=L)/(n*L)

ggplot(aes(x=x),data=data.frame(x=1:40))+
  stat_function(aes(color=" 2% (species/genus)"),fun=function(x) sharedkmerprop(30,0.01,x+1,10^6))+
  stat_function(aes(color=" 8% (species/genus)"),fun=function(x) sharedkmerprop(30,0.04,x+1,10^6))+
  stat_function(aes(color="10% (genus)"),fun=function(x) sharedkmerprop(30,0.05,x+1,10^6))+
  stat_function(aes(color="15% (family/genus)"),fun=function(x) sharedkmerprop(30,0.075,x+1,10^6))+
  stat_function(aes(color="20% (order/family)"),fun=function(x) sharedkmerprop(30,0.1,x+1,10^6))+
  stat_function(aes(color="30% (class/order)"),fun=function(x) sharedkmerprop(30,0.15,x+1,10^6))+
  stat_function(aes(color="40% (class/order)"),fun=function(x) sharedkmerprop(30,0.2,x+1,10^6))+
  stat_function(aes(color="40% (class/order)"),fun=function(x) 1)+
  theme_cowplot() +
  scale_x_continuous(trans = "log2")



pr = function(x,k,p,h=0)
  vapply(x,
         function(di) sum(choose(k,0:p)*di^(0:p)*(1-di)^(k-0:p)*(2*(1-h/k)^(0:p)-((1-h/k)^(0:p))^2)),
         c(1))

pr(0.15,30,4)

(pr(0.1,30,0)*vote(30,0)+
    pr(0.1,30,1)*vote(30,1)+
    pr(0.1,30,2)*vote(30,2)+
    pr(0.1,30,3)*vote(30,3)+
    pr(0.1,30,4)*vote(30,4)+
    pr(0.1,30,5)*vote(30,5))*120/10


ggplot()+stat_function(fun=function(x) pr(x,30,14,5))+
  stat_function(fun=function(x) pr(x,32,15,5),color="red")+
  coord_cartesian(xlim=c(0,.635))

vote = function(k,d) (1-d/k)^k

vote(30,4)*2

ggplot()+stat_function(fun=function(x) vote(30,x*5))+
  stat_function(fun=function(x) vote(32,x*5),color="red")+
  coord_cartesian(xlim=c(0,1))

require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

taxon_dist <- read_tsv("../data/closest_taxon_wrank.txt")

ggplot(taxon_dist) + aes(x=rank, y=dist, color=rank) +
  geom_boxplot() +
  geom_point()

changed = function(k,d) 1-(1-d)^k
atleasttwosame =  function(k,d,n=n) 1 - (dbinom(n-1,n,changed(k,d)) + dbinom(n,n,changed(k,d)))
sharedkmers = function(k,d,n=n,L=L) (n*L-n*L*atleasttwosame(k,d,n=n))/(n*L)


pr = function(x,k,p,h=0)
  vapply(x,
         function(di) sum(choose(k,0:p)*di^(0:p)*(1-di)^(k-0:p)*(2*(1-h/k)^(0:p)-((1-h/k)^(0:p))^2)),
         c(1))

pr(0.15,30,4)

(pr(0.1,30,0)*vote(30,0)+
    pr(0.1,30,1)*vote(30,1)+
    pr(0.1,30,2)*vote(30,2)+
    pr(0.1,30,3)*vote(30,3)+
    pr(0.1,30,4)*vote(30,4)+
    pr(0.1,30,5)*vote(30,5))*120/10


ggplot()+stat_function(fun=function(x) pr(x,30,14,5))+
  stat_function(fun=function(x) pr(x,32,15,5),color="red")+
  coord_cartesian(xlim=c(0,.635))

vote = function(k,d) (1-d/k)^k

vote(30,4)*2

ggplot()+stat_function(fun=function(x) vote(30,x*5))+
  stat_function(fun=function(x) vote(32,x*5),color="red")+
  coord_cartesian(xlim=c(0,1))