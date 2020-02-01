clear all

rng(1);

timeline = 10000000;
%general economy
a = 0.8;
f0r = 30;
f0 = f0r*(1-a);
fsig2 = 10;
sig = fsig2*(1-a^2);

%forestry
lf = 0.7; %lambda
af = 1;
ff0r = 25;
ff0 = ff0r*(1-lf*af);
ffsig2 = 2;
sigf = ffsig2*(1-lf^2*af^2);
factorf = lf*(1-af)*f0r/(1-af*lf);

%housing
lr = 0.7;
ar = 0;
fr0r = 25;
fr0 = fr0r*(1-lr*ar);
frsig2 = 2;
sigr = frsig2*(1-lr^2*ar^2);

%donation
lb = 0.7;
ab = 0.9;
fb0r = 25;
fb0 = fb0r*(1-lb*ab);
fbsig2 = 2;
sigb = fbsig2*(1-lb^2*ab^2);

% buy strategy stuff
code = [2];
al = 1;
be = 1;

godsimnum = 1;

period = 10;
lag = 0;
A = 1; % amplitude

burnin = 201; % burn in first few values of net return as they have not converged yet

cvalth = 0; % buying threshold for buystrat code 1

simtime = 9500000;%9500; % number of buying opportunities
fund = 0; % money saved
cumb = 0; % cummulative conservation value

bfn = 1; % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2

rho = 0.1; %0.1; % economic discount rate
del = 0.01; % ecological discount rate

efmu = 0; % mean of ef
efsig2 = 0; % var of ef
ermu = 0; % mean of er
ersig2 = 0; % var of er

ch = 0; % land change cost and option value

b_def = 10; % default b

param = [timeline,a,f0r,f0,fsig2,sig,lf,af,ff0r,ff0,ffsig2,sigf,factorf,lr,ar,fr0r,fr0,frsig2,sigr,...
lb,ab,fb0r,fb0,fbsig2,sigb,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];

receptacle = mainsim(param,code);
strcumb = receptacle{1};
fprintf("str     mean cumb\n");
stratstr = {'CVAL','Lc','Hc','Lff','Hff','Lfr','Hfr','LE','HE'};
for i = 1:length(strcumb)
  fprintf("%s       %.2f\n",stratstr{i},strcumb(i));
end