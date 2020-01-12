set.seed(5)
rm(list=ls())
pw= c(600,625) #pwot window

#general economy
t = seq(1,800)
a = 0.8
f0r = 30
f0 = f0r*(1-a)
sig = 10*(1-a^2)

#forestry
lf = 0.7 #lambda
af = 0.4
ff0r = 25
ff0 = ff0r*(1-lf*af)
sigf = 4*(1-lf^2*af^2)
factorf = lf*(1-af)*f0r/(1-af*lf)
#housing
lr = 0.7
ar = 0.1
fr0r = 25
fr0 = fr0r*(1-lr*ar)
sigr = 2*(1-lr^2*ar^2)

#donation
lb = 0.7
ab = 0.9
fb0r = 25
fb0 = fb0r*(1-lb*ab)
sigb = 4*(1-lb^2*ab^2)




for ( j in 1)
{
  f = rep(0,length(t))
  for ( i in 1:length(t))
  {
    if(i == 1)
    {
      f[i] = (f0 + rnorm(1))
    }
    else
    {
      f[i] = a*f[i-1] + (f0 + sig*rnorm(1))
    }
  }
  ff = rep(0,length(t))
  fr = rep(0,length(t))
  fb = rep(0,length(t))
  for (i in 1:length(t))
  {
    if( i == 1)
    {
      ff[i] = (ff0 + sigf*rnorm(1))
      fb[i] = (fb0 + sigb*rnorm(1))
      fr[i] = (fr0 + sigr*rnorm(1))
    }
    else
    {
      ff[i] = lf*(af*ff[i-1] + (1-af)*(f[i] - f0r)) + (ff0 + sigf*rnorm(1))
      fb[i] = lb*(ab*fb[i-1] + (1-ab)*(f[i] - f0r)) + (fb0 + sigb*rnorm(1))
      fr[i] = lr*(ar*fr[i-1] + (1-ar)*(f[i] - f0r)) + (fr0 + sigr*rnorm(1))
    }
  }
  
  if (j == 1)
  {
    plot(ff[pw[1]:pw[2]]~t[pw[1]:pw[2]],type='l',ylim=c(15,60))
    lines(fr[pw[1]:pw[2]]~t[pw[1]:pw[2]],col='blue')
    lines(fb[pw[1]:pw[2]]~t[pw[1]:pw[2]],col='green')
    lines(f[pw[1]:pw[2]]~t[pw[1]:pw[2]],col='red')
    legend('topright',legend=c('ff','fr','fb','f'),col=c('black','blue','green','red'),lty=1)
  }
  else
  {
    lines(ff[pw[1]:pw[2]]~t[pw[1]:pw[2]])
    lines(f[pw[1]:pw[2]]~t[pw[1]:pw[2]],col='red')
  }
}

sprintf('analytic ff mean = %f',(lf*(1-af)*f0r)/(1-af*lf) + ff0r) #analytic ff mean
sprintf('simulation ff mean = %f',mean(ff[200:800])) #simulation mean

sprintf('analytic f mean = %f',f0r) #analytic ff mean
sprintf('simulation f mean = %f',mean(f[200:800])) #simulation mean

sprintf("cor(ff,f)= %.2f",cor(ff,f))
sprintf("cor(fr,f)= %.2f",cor(fr,f))
sprintf("cor(fb,f)= %.2f",cor(fb,f))
