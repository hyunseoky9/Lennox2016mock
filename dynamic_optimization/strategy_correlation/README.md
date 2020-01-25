looking at correlation between strategy rewards and autocorrelation coefficients, and cross correlation btw strategies' rewards with a set coefficients.

5scenarios:
From default and simtime=50,
1. When cor(f,fr) ~ 1, ff cahgnes from low cor(f,ff) to high cor(f,ff)
param: lr=0.7, ar=0; lf=0.8, af=1~0
2. Opposite of 1. (ff is now fr)
param: lf=0.7, af=0; lr=0.8, ar=1~0
3. fb0r goes from low to high 
param: fb0r = fr0r ~ E(C)
4. cor(fb,f) goes from low to high, when fb0r ~ fr0r
param: fb0r = fr0r, ab = 1~0
5. same as 4. but when fb0r ~ E(C)
param: lb=0.8, ab=1~0

