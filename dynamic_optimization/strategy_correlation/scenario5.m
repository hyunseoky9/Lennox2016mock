clear all

%rng(1);
timeline = 800;
%general economy
a = 0.8;
x0r = 30;
x0 = x0r*(1-a);
xsig2 = 10;
sig = xsig2*(1-a^2);

code = 1:9;

ab = linspace(1,0,11);
M = zeros(length(ab),1+length(code));
M(:,1) = ab;
for i = 1:length(ab)
	fprintf("ab=%.2f\n",ab(i));
	%forestry
	lf = 0.7; %lambda
	af = 0.4; %linspace(1,0,101);
	xf0r = 25;
	xf0 = xf0r*(1-lf*af);
	xfsig2 = 2;
	sigf = xfsig2*(1-lf^2*af^2);
	factorf = lf*(1-af)*x0r/(1-af*lf);

	%housing
	lr = 0.7;
	ar = 0.4;
	xr0r = 25;
	xr0 = xr0r*(1-lr*ar);
	xrsig2 = 2;
	sigr = xrsig2*(1-lr^2*ar^2);

	%donation
	lb = 0.7;
	xb0r = 275;
	xb0 = xb0r*(1-lb*ab(i));
	xbsig2 = 2;
	sigb = xbsig2*(1-lb^2*ab(i)^2);

	% buy strategy stuff
	al = 1;
	be = 1;

	godsimnum = 100;

	period = 10;
	lag = 0;
	A = 1; % amplitude

	burnin = 201; % burn in first few values of net return as they have not converged yet

	cvalth = 0; % buying threshold for buystrat code 1

	simtime = 50; % number of buying opportunities
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



	param = [timeline,a,x0r,x0,xsig2,sig,lf,af,xf0r,xf0,xfsig2,sigf,factorf,lr,ar,xr0r,xr0,xrsig2,sigr,...
	lb,ab(i),xb0r,xb0,xbsig2,sigb,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
	cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];
	receptacle = mainsim(param,code);
	strcumb = receptacle{1};
	M(i,2:end) = strcumb;
end

fprintf("cumb\n");
fprintf("    af        CVAL       LC        Hc       Lxf       Hxf      Lxr      Hxr        LE        HE\n");
display(M);

fprintf("corr coeff\n");
fprintf("    af        CVAL       LC        Hc       Lxf       Hxf      Lxr      Hxr        LE        HE\n");
display(corrcoef(M));
%fprintf("str     mean cumb\n");
%stratstr = {'CVAL','Lc','Hc','Lxf','Hxf','Lxr','Hxr','LE','HE'};
%for i = 1:length(strcumb)
%  fprintf("%s       %.2f\n",stratstr{i},strcumb(i));
%endn