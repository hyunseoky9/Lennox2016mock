clear all
t = 100;
a = 0.8;
x0r = 30;
xsig2 = 10;
x0 = x0r*(1-a);
sig = xsig2*(1-a^2);

%forestry
lf = 0.9; %lambda
af = 1;
xf0r = 50;
xfsig2 = 2;
xf0 = xf0r*(1-lf*af);
sigf = xfsig2*(1-lf^2*af^2);
e = -3;
efj = e;

%housing
lr = 0.9;
ar = 1;
xr0r = 50;
xrsig2 = 2;
xr0 = xr0r*(1-lr*ar);
sigr = xrsig2*(1-lr^2*ar^2);
erj = -e;

Ex = zeros([1 t]);
Exf = zeros([1 t]);
Exr = zeros([1 t]);

h = 10;
xst = 33; %x start value
xfst = xf0r + h; %xf start value
xrst = xr0r - h; %xr start value

val = 0;
valf = 0;
valr = 0;
for i = 1:t
  val = a^i*xst;
  valf = lf^i*af^i*xfst + efj;
  valr = lr^i*ar^i*xrst + erj;
  for j = 0:i-1
    val = val + a^j*x0;
    valf = valf + lf^j*af^j*xf0;
    valr = valr + lr^j*ar^j*xr0;
  end
  Ex(i) = val;
  Exf(i) = valf;
  Exr(i) = valr;
end
xfst = xfst+efj;
xrst = xrst+erj;
Exf = [xfst Exf];
Exr = [xrst Exr];


% plotting 
subplot(2,1,1);
plot(0:t,Exf,'-');
hold on   
plot(0:t,Exr,'--');
legend('E(xfj)','E(xrj)');
xlabel('t') 
hold off 


% getting PV
del = 0.99;
d = del.^(1:t);
PVxf = zeros([1 t]);
PVxf(t+1) = Exf(t+1)*del^t;
PVxr(t+1) = Exr(t+1)*del^t;
for i = t:-1:1  
  PVxf(i) = PVxf(i+1) + Exf(i)*del^(i-1);
  PVxr(i) = PVxr(i+1) + Exr(i)*del^(i-1);
end

% plotting PV
%subplot(3,1,2);
%plot(0:t,PVxf);
%hold on
%plot(0:t,PVxr);
%hold off
%legend('PVxfj','PVxrj');

pvdiff = PVxf-PVxr;
subplot(2,1,2);
plot(0:t,pvdiff);
xlabel('tj')
legend("PVxfj-PVxrj");

[val,ind] = min(pvdiff);
[val,ind2] = find(Exf-Exr < 0);
fprintf("min(PVxfj-PVxfj) occur when t=%d\n",ind-1);
fprintf("when Exfj-Exrj becomes < 0  t=%d\n",ind2(1)-1);