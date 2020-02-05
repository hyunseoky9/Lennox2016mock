looking at correlation between strategy rewards and autocorrelation coefficients, and cross correlation btw strategies' rewards with a set coefficients.

5scenarios:
From default parameter (look at default_parameter_values.txt) values and simtime=50,
1. When cor(x,xr) ~ 1, xf cahgnes from low cor(x,xf) to high cor(x,xf)
param: lr=0.7, ar=0; lf=0.8, af=1~0
2. Opposite of 1. (xf is now xr)
param: lf=0.7, af=0; lr=0.8, ar=1~0
3. xb0r goes from low to high 
param: xb0r = xr0r ~ E(C)
4. cor(xb,x) goes from low to high, when xb0r ~ xr0r
param: xb0r = xr0r, ab = 1~0
5. same as 4. but when xb0r ~ E(C)
param: lb=0.8, ab=1~0

