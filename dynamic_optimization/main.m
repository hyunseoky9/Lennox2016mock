%clear all

timeline = 10000; %10000;
%general economy
a = 0.8;
x0r = 30;
x0 = x0r*(1-a);
xsig2 = 10;
sig = xsig2*(1-a^2);


%forestry
lf = 0.7; %lambda
af = 0;
xf0r = 25;
xf0 = xf0r*(1-lf*af);
xfsig2 = 2;
sigf = xfsig2*(1-lf^2*af^2);
factorf = lf*(1-af)*x0r/(1-af*lf);

%housing
lr = 0.7;
ar = 1;
xr0r = 25;
xr0 = xr0r*(1-lr*ar);
xrsig2 = 2;
sigr = xrsig2*(1-lr^2*ar^2);

%donation
lb = 0.7;
ab = 0.9;
xb0r = 25;
xb0 = xb0r*(1-lb*ab);
xbsig2 = 2;
sigb = xbsig2*(1-lb^2*ab^2);

% buy strategy stuff
code = [2,3];
al = 1;
be = 1;

godsimnum = 1;

period = 10;
lag = 0;
A = 1; % amplitude

burnin = 201; % burn in first few values of net return as they have not converged yet

cvalth = 0; % buying threshold for buystrat code 1

simtime = 9500; % number of buying opportunities
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

parambundle = {timeline,a,x0r,xsig2,lf,af,xf0r,xfsig2,lr,ar,xr0r,xrsig2,...
lb,ab,xb0r,xbsig2,code,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def};
paramset = paramsetmaker(parambundle);

stratstr = {'CVAL','Lc','Hc','Lxf','Hxf','Lxr','Hxr','LE','HE'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter changes
%for i = 1:2
%	if i == 1
%		af = 1;
%		ar = 0;
%		param = [timeline,a,x0r,x0,xsig2,sig,lf,af,xf0r,xf0,xfsig2,sigf,factorf,lr,ar,xr0r,xr0,xrsig2,sigr,...
%		lb,ab,xb0r,xb0,xbsig2,sigb,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
%		cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];
%		receptacle = mainsim(param,code);
%		figure(i)
%		bleh = plotting(receptacle,param);
%	else 
%		af = 0;
%		ar = 1;
%		param = [timeline,a,x0r,x0,xsig2,sig,lf,af,xf0r,xf0,xfsig2,sigf,factorf,lr,ar,xr0r,xr0,xrsig2,sigr,...
%		lb,ab,xb0r,xb0,xbsig2,sigb,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
%		cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];
%		receptacle2 = mainsim(param,code);
%		figure(i)
%		bleh = plotting(receptacle2,param);
%	end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(paramset);
  param = paramset{i};
  receptacle = mainsim(param);
  strcumb = receptacle{1};
end

param = [timeline,a,x0r,xsig2,lf,af,xf0r,xfsig2,lr,ar,xr0r,xrsig2,...
lb,ab,xb0r,xbsig2,code(1),al(1),be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];
receptacle = mainsim(param);
bleh = plotting(receptacle,param);

strcumb = receptacle{1};
t = receptacle{2};
x = receptacle{3};
xf = receptacle{4};
xr = receptacle{5};
xb = receptacle{6};
tjs = receptacle{7};
C = receptacle{8};
fundt = receptacle{9};
ben = receptacle{10};
buy = receptacle{11};

%if mod(god,100) == 0
%	fprintf('god=%d\n',god);
%end
fprintf("ar=%.2f, af=%.2f\n",ar,af);
fprintf("str     mean cumb\n");
for i = 1:length(strcumb)
  fprintf("%s       %.2f\n",stratstr{code(i)},strcumb(i));
end