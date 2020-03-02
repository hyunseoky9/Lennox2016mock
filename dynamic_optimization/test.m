<<<<<<< HEAD
a = 0.8;
x0r = 30;
xsig2 = 10;
x0 = x0r*(1-a);
sig = xsig2*(1-a^2);

%forestry
lf = 0.9; %lambda
af = 1;
xf0r = 25;
xfsig2 = 2;
xf0 = xf0r*(1-lf*af);
sigf = xfsig2*(1-lf^2*af^2);

Ex = zeros([1 100]);
Exf = zeros([1 100]);

xst = 33; 
xfst = 1;

val = 0;
valf = 0;
for i = 1:100
  val = a^i*xst;
  valf = lf^i*af^i*xfst;
  for j = 0:i-1
    val = val + a^j*x0;
    valf = valf + lf^j*af^j*xf0;
  Ex(i) = val;
  Exf(i) = valf;
  end
end

plot(1:100,Exf);
=======
efmu = 0;
efsig = 1;
mu = [efmu 10];
sigma = [1 0.8; 0.8 1];
R = mvnrnd(mu,sigma,1);
e_fj = R(1); % indiv forest var
b = R(2); % for benefit
R
>>>>>>> parent of 9232945... its been awhile...
