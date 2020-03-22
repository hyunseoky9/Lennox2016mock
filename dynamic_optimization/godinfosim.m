function [receptacle] = godinfosim(param)
%mxtjcor = 0;
%mxCcor = 0;
%mxBcor = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(1);

stratstr = {'CVAL','Lc','Hc','Lxf','Hxf','Lxr','Hxr','LE','HE','HB','LB','Hx','Lx'};

pind = 1; % parameter index
timeline = param(pind);
pind = pind + 1;
t = linspace(1,timeline,timeline);

%general economy
a = param(pind);
pind = pind + 1;
x0r = param(pind);
pind = pind + 1;
x0 = x0r*(1-a);
xsig2 = param(pind);
pind = pind + 1;
sig = xsig2*(1-a^2);

%forestry
lf = param(pind); %lambda
pind = pind + 1;
af = param(pind);
pind = pind + 1;
xf0r = param(pind);
pind = pind + 1;
xf0 = xf0r*(1-lf*af);
xfsig2 = param(pind);
pind = pind + 1;
sigf = xfsig2*(1-lf^2*af^2);
factorf = lf*(1-af)*x0r/(1-af*lf);

%housing
lr = param(pind);
pind = pind + 1;
ar = param(pind);
pind = pind + 1;
xr0r = param(pind);
pind = pind + 1;
xr0 = xr0r*(1-lr*ar);
xrsig2 = param(pind);
pind = pind + 1;
sigr = xrsig2*(1-lr^2*ar^2);

%donation
lb = param(pind);
pind = pind + 1;
ab = param(pind);
pind = pind + 1;
xb0r = param(pind);
pind = pind + 1;
xb0 = xb0r*(1-lb*ab);
xbsig2 = param(pind);
pind = pind + 1;
sigb = xbsig2*(1-lb^2*ab^2);

% buy strategy stuff
al = param(pind);
pind = pind + 1;
be = param(pind);
pind = pind + 1;

godsimnum = param(pind);
pind = pind + 1;
godsim = godsimnum;

% params for deterministic periodic function for xf,xr,xb
period = param(pind);
pind = pind + 1;
lag = param(pind);
pind = pind + 1;
A = param(pind); % amplitude
pind = pind + 1;

burnin = param(pind); % burn in first few values of net return as they have not converged yet
pind = pind + 1;

cvalth = param(pind); % buying threshold for buystrat code 1
pind = pind + 1;

simtime = param(pind); % number of buying opportunities
pind = pind + 1;
fund = param(pind); % money saved
pind = pind + 1;
fundt = zeros(1,simtime); % array for fund over time
cumb = param(pind); % cummulative conservation value
pind = pind + 1;

bfn = param(pind); % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2
pind = pind + 1;
rho = param(pind); %0.1; % economic discount rate
pind = pind + 1;
del = param(pind); % ecological discount rate
pind = pind + 1;

efmu = param(pind); % mean of ef
pind = pind + 1;
efsig2 = param(pind); % var of ef
pind = pind + 1;
ermu = param(pind); % mean of er
pind = pind + 1;
ersig2 = param(pind); % var of er
pind = pind + 1;

ch = param(pind); % land change cost and option value
pind = pind + 1;

b_def = param(pind); % default b
pind = pind + 1;
t_jmethod = 1;
pind = pind + 1;
intrate = 0;
pind = pind + 1;

code = param(pind);
pind = pind + 1;


Exf = xf0r;
Vxf = xfsig2 + lf^2*(1-af)*xsig2/(1-af^2*lf^2);
Exr = xr0r;
Vxr = xrsig2 + lr^2*(1-ar)*xsig2/(1-ar^2*lr^2);
Lxf = Exf; %norminv(1/3,Exf,Vxf^(1/2)); % threshold for low xf
Hxf = Exf; %norminv(2/3,Exf,Vxf^(1/2)); % threshold for high xf
Lxr = Exr; %norminv(1/3,Exr,Vxr^(1/2)); % threshold for low xr
Hxr = Exr; %norminv(2/3,Exr,Vxr^(1/2)); % threshold for high xr
Lx = x0r; %norminv(1/3,x0r,xsig2^(1/2)); % threshold for low x
Hx = x0r; %norminv(2/3,x0r,xsig2^(1/2)); % threshold for high x

%tjsrecep = zeros(1,param(pind+));
strcumb = 0; % strategy's cumulative benefit
for god = 1:godsim
  x = zeros(0,length(t));
  xf = zeros(1,length(t));
  xr = zeros(1,length(t));
  xb = zeros(1,length(t));
  regsin = 0; % 0=autoregression 1=sinusoidal (for xf,xr,xb)
  if regsin == 0   
    nr = normrnd(0,1,[1 length(t)]);
    x(1) = (x0 + normrnd(0,1));
    for  i = 2:length(t)
      x(i) = a*x(i-1) + (x0 + sig*nr(i)); %normrnd(0,1));
    end

    nr = normrnd(0,1,[1 length(t)]);
    xf(1) = (xf0 + sigf*normrnd(0,1));
    xb(1) = (xb0 + sigb*normrnd(0,1));
    xr(1) = (xr0 + sigr*normrnd(0,1));
    for i = 2:length(t)
      xf(i) = lf*(af*xf(i-1) + (1-af)*(x(i) - x0r)) + (xf0 + sigf*nr(i));%normrnd(0,1));
    end
    nr = normrnd(0,1,[1 length(t)]);
    for i = 2:length(t)
      xb(i) = lb*(ab*xb(i-1) + (1-ab)*(x(i) - x0r)) + (xb0 + sigb*nr(i));%normrnd(0,1));
    end
    nr = normrnd(0,1,[1 length(t)]);
    for i = 2:length(t)
      xr(i) = lr*(ar*xr(i-1) + (1-ar)*(x(i) - x0r)) + (xr0 + sigr*nr(i));%normrnd(0,1));
    end
  else
    x = x0r + A*cos(2*pi*(1:length(t))/period-pi);
    xf = xf0r + A*cos(2*pi*(1:length(t))/period);
    xr = xr0r + A*cos(2*pi*(1:length(t))/period-pi);
    xb = xb0r + A*cos(2*pi*(1:length(t))/period-pi);
  end



  %fprintf('\n');
  %fprintf('analytic xf mean = %f\n', xf0r) %analytic xf mean
  %fprintf('simulation xf mean = %f\n',mean(xf(200:800))) %simulation mean

  %fprintf('analytic x mean = %f\n',x0r) %analytic xf mean
  %fprintf('simulation x mean = %f\n',mean(x(200:800))) %simulation mean
  %fprintf('\n');
  %fprintf("correlations btw net returns\n");
  %cor = corrcoef(xf,x);
  %fprintf("cor(xf,x)= %.2f\n",cor(1,2))
  %cor = corrcoef(xr,x);
  %fprintf("cor(xr,x)= %.2f\n",cor(1,2))
  %cor = corrcoef(xb,x);
  %fprintf("cor(xb,x)= %.2f\n",cor(1,2))
  %fprintf("\n");

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % simulation of tj, b, c, buying, etc.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  x = x(burnin:end);
  xf = xf(burnin:end);
  xb = xb(burnin:end);
  xr = xr(burnin:end);
  
  edisc = repelem(1+rho,length(x)).^(0:(length(x)-1)); %economic discount rate.
  ecodisc = repelem(1+del,length(x)).^(0:(length(x)-1)); %ecological discount rate

  C = zeros(1,simtime);
  ben = zeros(1,simtime);
  tjs = zeros(1,simtime);
  buy = [];
  %mcval = [Inf 0];
  %cvalm = 0;
  %% tj and cost 
  for i = 1:simtime
    %fprintf('time step=%d\n',i);
    
    nxf = xf(i:end);
    nxb = xb(i:end);
    nxr = xr(i:end);
    %fprintf('length of xb=%d\n',size(xb,2));
    %fprintf('this yrs donation = %.2f\n',xb(1));
    % getting e_fj,e_rj, and b (if bfn=2)
   

    e_rj = normrnd(ermu,ersig2); % indiv dev var

    if bfn == 1
      e_fj = normrnd(efmu,efsig2);
      b = b_def;
    elseif bfn >= 2
      b_mu = b_def;
      mu = [efmu b_mu];
      sigma = [efsig2 0.8; 0.8 1];
      R = mvnrnd(mu,sigma,1);
      e_fj = R(1); % indiv forest var
      b = R(2); % for benefit
    end
    if bfn == 3
      b = b^(1/2);
    end
    x_fj = nxf + e_fj; 
    x_fj(x_fj<0) = 0; % j's xf
    x_rj = nxr + e_rj; 
    x_rj(x_rj<0) = 0; % j's xr


    % getting clearing time t_j
    t_j = 1;
    nedisc = edisc(1:end-i+1); % new economic discount rate for this sim.
    if t_jmethod
      PVdiff = sum((x_rj(t_j:end)-x_fj(t_j:end))./nedisc(1:(end-t_j+1))) - ch;
      PVdiffarray = zeros(1,(length(x_rj)-t_j+1));
      PVdiffarray(1) = PVdiff;
      for oi = 2:length(PVdiffarray)
        PVdiffarray(oi) = PVdiffarray(oi-1) - (x_rj(oi-1)-x_fj(oi-1))/nedisc(oi-1);
      end

      [val, t_j] = max(PVdiffarray);
      if t_j > length(x_fj)
        t_j = length(x_fj);
      end
      if val < 0
        t_j = length(x_fj);
      end
    else
      PVdiff = sum((x_rj(t_j:end)-x_fj(t_j:end))./nedisc(1:(end-t_j+1))) - ch;
      while PVdiff <= 0 && t_j <= length(x_fj)
        PVdiff = PVdiff - (x_rj(t_j)-x_fj(t_j))/nedisc(t_j);
        t_j = t_j + 1;
      end
    end

    tjs(i) = t_j;
    %fprintf("tj=%d\n",t_j);
    if t_j > 1
        %fprintf('tj was bigger than 1 on step %d; t_j=%d\n',i,t_j);
    end
    %cind = 1;
    %t_j = cind; %% test
    % cost
    if t_j == 1
      c = sum(x_rj(t_j:end)./nedisc(t_j:end));
      %c = sum(f_fj(t_j:end)./nedisc(t_j:end));
    else
      %c = sum(f_rj(1:end)./nedisc(1:end));
      %c = sum(f_fj(1:end)./nedisc(1:end));
      %foo = 2;
      %c = sum(f_fj(1:(foo-1))./nedisc(1:(foo-1))) + sum(f_rj(foo:end)./nedisc(foo:end));
      %c = sum(f_rj(1:(foo-1))./nedisc(1:(foo-1))) + sum(f_fj(foo:end)./nedisc(foo:end));
      c = sum(x_fj(1:(t_j-1))./nedisc(1:(t_j-1))) + sum(x_rj(t_j:end)./nedisc(t_j:end));
    end

    %C(i) = c;

    %% only getting C with certain t_j
    %cind = 1;
    %Ctitle = sprintf('tj = %d',cind);
    %if t_j == cind
    C(i) = c;
    %end

    %v = f_rj(t_j:end)./edisc(t_j:end);
    %disp(v(1:15))
  %% benefit and buying
    %benefit
    necodisc = ecodisc(1:end-i+1); %ecological discount rate
    B = sum(b./necodisc(t_j:end)); % conservation value (incl' discount and threat component)
    ben(i) = B;
    %v = b ./ecodisc(t_j:end);
    %disp(v(1:15));
  end
  %tjstemp = tjs;
  %tjstemp(tjstemp>1) = 0;
  %tjsrecep = tjsrecep + tjstemp;
  %fprintf("min(C)=%.2f, max(C)=%.2f\n",min(C),max(C));
  E = ben./C; % roi ratio
  
  %% evaluating and buying process
  nenough = 0;
  ncrinotmet = 0;
  mc = mean(C);
  mB = mean(ben);
  me = mean(E);
  for i = 1:simtime
    fund = fund + xb(i); % add this yr's fund to the account
    %[cumb,fund,buy] = buystrat(buy,code(8),cumb,fund,xb(i),ben(i),C(i),E(i),xf(i),xr(i),al,be,cvalth,Lc,Hc,Lxf,Hxf,Lxr,Hxr,Lx,Hx,LE,HE);
    [cumb,fund,buy,notenough,crinotmet] = buystrat(buy,code,cumb,fund,xb(i),ben(i),C(i),E(i),x(i),xf(i),xr(i),al,be,cvalth,mc,Lxf,Hxf,Lxr,Hxr,Lx,Hx,me,mB);
    fund = fund*(1+intrate);
    fundt(i) = fund;
    nenough = nenough + notenough;
    ncrinotmet = ncrinotmet + crinotmet;
  end
  strcumb = strcumb + cumb;
  %fprintf('cumb=%.2f\n',cumb);
  %fprintf('bought %d times. criteria not met %d times. fund not enough %d times. with code %s\n',sum(buy),ncrinotmet,nenough,stratstr{code});
  cumb = 0;
  %% calculate correlation btw f and C,ben, tjs.
  %fprintf("corelations between f and C,B,tj\n");
  %cor = corrcoef(x(1:simtime),C);
  %mxCcor = mxCcor + cor(1,2);
  %fprintf('cor(x,C)=%.2f\n',cor(1,2));
  %cor = corrcoef(x(1:simtime),ben);
  %mxBcor = mxBcor + cor(1,2);
  %fprintf('cor((x(1:simtime),tjs);
  %cor = corrcoef(x(1:simtime),tjs);
  %mxtjcor = mxtjcor + cor(1,2);
  %fprintf('cor(x,tj)=%.2f\n',cor(1,2));
  %fprintf('\n');

  %% calculating tj ratio
  %y = zeros(2,length(unique(tjs)));
  %temp = unique(tjs);
  %y(1,:) = unique(tjs);
  %for i = 1:length(unique(tjs))   
  %  y(2,i) = sum(tjs==temp(i));
  %end
  %fprintf("tj,ratio\n");
  %for i = 1:length(unique(tjs))
  %  fprintf('%d,%.2f\n',y(1,i),y(2,i)/length(tjs));
  %end

  %fprintf('cumb = .%2f\n',cumb);
  %fprintf('mincval = %.2f, maxcval = %.2f\n',mcval(1),mcval(2));
  %fprintf("mean cval = %.2f\n",cvalm/simtime);
  %fprintf('buy=%d \n',find(buy > 0));

  if mod(god,100) == 0
    if god == 100
      nchar = fprintf('god=%d',god);
    elseif god >= 100
      fprintf(repmat('\b', 1, nchar));
      fprintf('god=%d',god);
    end
  end
end
fprintf('\n');

%tjsrecep = tjsrecep/(godsim/length(code));
%plot(1:length(tjsrecep),tjsrecep);
%xlabel('tjsrecep');

strcumb = strcumb/godsim; % average over sim rep
receptacle = {strcumb,t,x,xf,xr,xb,tjs,C,fundt,ben,buy,nenough,ncrinotmet};
%fprintf("str     mean cumb\n");
%for i = 1:length(strcumb)
%  fprintf("%s       %.2f\n",stratstr{i},strcumb(i));
%end


% fprintf('mean of cor(x,C)=%.2f\n',mxCcor/godsim);
% fprintf('mean of cor(x,B)=%.2f\n',mxBcor/godsim);
% fprintf('mean of cor(x,tj)=%.2f\n',mxtjcor/godsim);
%fprintf("mean tj=%.2f\n",sum(tjs)/simtime);
%fprintf("pr(tj=1) = %.2f\n", length(tjs(tjs==1))/length(tjs));
%fprintf("pr(tj=2) = %.2f\n", length(tjs(tjs==2))/length(tjs));
%fprintf("pr(tj=3) = %.2f\n", length(tjs(tjs==3))/length(tjs));