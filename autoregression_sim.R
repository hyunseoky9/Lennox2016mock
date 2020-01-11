t = seq(1,800)
a = 0.8
f0 = 10*(1-a)
sig = (1-a^2)


lf = 0.7 #lambda
af = 0.8
ff0 = 5*(1-lf*af)
sigf = (1-lf^2*af^2)

for ( j in 1:20)
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
  for (i in 1:length(t))
  {
    if( i == 1)
    {
      ff[i] = (ff0 + sigf*rnorm(1))
    }
    else
    {
      ff[i] = lf*(af*ff[i-1] + (1-af)*f[i-1]) + (ff0 + sigf*rnorm(1))
    }
  }
  
  if (j == 1)
  {
    plot(ff~t,type='l',ylim=c(0,15))
    lines(f~t,col='red')
  }
  else
  {
    lines(ff~t)
    lines(f~t,col='red')
  }
}

sprintf('analytic ff mean = %f',(lf*(1-af)*10)/(1-af*lf) + 5) #analytic ff mean
sprintf('simulation ff mean = %f',mean(ff[200:800])) #simulation mean

sprintf('analytic f mean = %f',10) #analytic ff mean
sprintf('simulation f mean = %f',mean(f[200:800])) #simulation mean
