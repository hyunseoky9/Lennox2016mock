function [receptacle] = mainsim(param,code)
mftjcor = 0;
mfCcor = 0;
mfBcor = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pind = 1; % parameter index
%rng(1);
timeline = param(pind);
t = linspace(1,timeline,timeline);

%general economy
a = param(pind+1);
f0r = param(pind+2);
f0 = param(pind+3);
fsig2 = param(pind+4);
sig = param(pind+5);

%forestry
lf = param(pind+6); %lambda
af = param(pind+7);
ff0r = param(pind+8);
ff0 = param(pind+9);
ffsig2 = param(pind+10);
sigf = param(pind+11);
factorf = param(pind+12);

%housing
lr = param(pind+13);
ar = param(pind+14);
fr0r = param(pind+15);
fr0 = param(pind+16);
frsig2 = param(pind+17);
sigr = param(pind+18);

%donation
lb = param(pind+19);
ab = param(pind+20);
fb0r = param(pind+21);
fb0 = param(pind+22);
fbsig2 = param(pind+23);
sigb = param(pind+24);

% buy strategy stuff
code = code;
al = param(pind+25);
be = param(pind+26);

godsimnum = param(pind+27);
godsim = length(code)*godsimnum;

% params for deterministic periodic function for ff,fr,fb
period = param(pind+28);
lag = param(pind+29);
A = param(pind+30); % amplitude

burnin = param(pind+31); % burn in first few values of net return as they have not converged yet

cvalth = param(pind+32); % buying threshold for buystrat code 1

simtime = param(pind+33); % number of buying opportunities
fund = param(pind+34); % money saved
fundt = zeros(1,simtime); % array for fund over time
cumb = param(pind+35); % cummulative conservation value

bfn = param(pind+36); % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2
rho = param(pind+37); %0.1; % economic discount rate
del = param(pind+38); % ecological discount rate

efmu = param(pind+39); % mean of ef
efsig2 = param(pind+40); % var of ef
ermu = param(pind+41); % mean of er
ersig2 = param(pind+42); % var of er

ch = param(pind+43); % land change cost and option value

b_def = param(pind+44); % default b


%tjsrecep = zeros(1,param(pind+));
strcumb = zeros(1,length(code)); % strategy's cumulative benefit
for god = 1:godsim
  f = zeros(0,length(t));
  ff = zeros(1,length(t));
  fr = zeros(1,length(t));
  fb = zeros(1,length(t));
  regsin = 0; % 0=autoregression 1=sinusoidal (for ff,fr,fb)
  if regsin == 0   
    f(1) = (f0 + normrnd(0,1));
    for  i = 2:length(t)
      f(i) = a*f(i-1) + (f0 + sig*normrnd(0,1));
    end

    ff(1) = (ff0 + sigf*normrnd(0,1));
    fb(1) = (fb0 + sigb*normrnd(0,1));
    f(1) = (fr0 + sigr*normrnd(0,1));
    for i = 2:length(t)
      ff(i) = lf*(af*ff(i-1) + (1-af)*(f(i) - f0r)) + (ff0 + sigf*normrnd(0,1));
      fb(i) = lb*(ab*fb(i-1) + (1-ab)*(f(i) - f0r)) + (fb0 + sigb*normrnd(0,1));
      fr(i) = lr*(ar*fr(i-1) + (1-ar)*(f(i) - f0r)) + (fr0 + sigr*normrnd(0,1));
    end
  else
    f = f0r + A*cos(2*pi*(1:length(t))/period-pi);
    ff = ff0r + A*cos(2*pi*(1:length(t))/period);
    fr = fr0r + A*cos(2*pi*(1:length(t))/period-pi);
    fb = fb0r + A*cos(2*pi*(1:length(t))/period-pi);
  end



  %fprintf('\n');
  %fprintf('analytic ff mean = %f\n', ff0r) %analytic ff mean
  %fprintf('simulation ff mean = %f\n',mean(ff(200:800))) %simulation mean

  %fprintf('analytic f mean = %f\n',f0r) %analytic ff mean
  %fprintf('simulation f mean = %f\n',mean(f(200:800))) %simulation mean
  %fprintf('\n');
  %fprintf("correlations btw net returns\n");
  %cor = corrcoef(ff,f);
  %fprintf("cor(ff,f)= %.2f\n",cor(1,2))
  %cor = corrcoef(fr,f);
  %fprintf("cor(fr,f)= %.2f\n",cor(1,2))
  %cor = corrcoef(fb,f);
  %fprintf("cor(fb,f)= %.2f\n",cor(1,2))
  %fprintf("\n");

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % simulation of tj, b, c, buying, etc.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  f = f(burnin:end);
  ff = ff(burnin:end);
  fb = fb(burnin:end);
  fr = fr(burnin:end);

  Eff = ff0r;
  Vff = ffsig2 + lf^2*(1-af)*fsig2/(1-af^2*lf^2);
  Efr = fr0r;
  Vfr = frsig2 + lr^2*(1-ar)*fsig2/(1-ar^2*lr^2);
  Lff = norminv(1/3,Eff,Vff^(1/2)); % threshold for low ff
  Hff = norminv(2/3,Eff,Vff^(1/2)); % threshold for high ff
  Lfr = norminv(1/3,Efr,Vfr^(1/2)); % threshold for low fr
  Hfr = norminv(2/3,Efr,Vfr^(1/2)); % threshold for high fr
  Lf = norminv(1/3,f0r,fsig2^(1/2)); % threshold for low f
  Hf = norminv(2/3,f0r,fsig2^(1/2)); % threshold for high f
  
  edisc = repelem(1+rho,length(f)).^(0:(length(f)-1)); %economic discount rate.
  ecodisc = repelem(1+del,length(f)).^(0:(length(f)-1)); %ecological discount rate

  C = zeros(1,simtime);
  ben = zeros(1,simtime);
  tjs = zeros(1,simtime);
  buy = [];
  mcval = [Inf 0];
  cvalm = 0;
  %% tj and cost 
  for i = 1:simtime
    %fprintf('time step=%d\n',i);
    
    nf = f(i:end);
    nff = ff(i:end);
    nfb = fb(i:end);
    nfr = fr(i:end);
    %fprintf('length of fb=%d\n',size(fb,2));
    %fprintf('this yrs donation = %.2f\n',fb(1));
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
    f_fj = nff + e_fj; 
    f_fj(f_fj<0) = 0; % j's ff
    f_rj = nfr + e_rj; 
    f_rj(f_rj<0) = 0; % j's fr


    % getting clearing time t_j
    t_jmethod = 0; % 0=earliest time its profitable, 1=time when its most profitable
    t_j = 1;
    if t_jmethod
      nedisc = edisc(1:end-i+1); % new economic discount rate for this sim.
      tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./nedisc(1:(end-t_j+1))) - ch;
      optimarray = zeros(1,(length(f_rj)-t_j+1));
      optimarray(1) = tjoptim;
      for oi = 2:length(optimarray)
        optimarray(oi) = optimarray(oi-1) - (f_rj(oi)-f_fj(oi))/nedisc(oi);
      end

      [val, t_j] = max(optimarray);
      if t_j > length(f_fj)
        t_j = length(f_fj);
      end
      if val < 0
        t_j = length(f_fj);
      end
    else
      nedisc = edisc(1:end-i+1); % new economic discount rate for this sim.
      tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./nedisc(1:(end-t_j+1))) - ch;
      while tjoptim <= 0 && t_j <= length(f_fj)
        tjoptim = tjoptim - (f_rj(t_j)-f_fj(t_j))/nedisc(t_j);
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
      c = sum(f_rj(t_j:end)./nedisc(t_j:end));
      %c = sum(f_fj(t_j:end)./nedisc(t_j:end));
    else
      %c = sum(f_rj(1:end)./nedisc(1:end));
      %c = sum(f_fj(1:end)./nedisc(1:end));
      %foo = 2;
      %c = sum(f_fj(1:(foo-1))./nedisc(1:(foo-1))) + sum(f_rj(foo:end)./nedisc(foo:end));
      %c = sum(f_rj(1:(foo-1))./nedisc(1:(foo-1))) + sum(f_fj(foo:end)./nedisc(foo:end));
      c = sum(f_fj(1:(t_j-1))./nedisc(1:(t_j-1))) + sum(f_rj(t_j:end)./nedisc(t_j:end));
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
  Eint = max(E) - min(E);
  LE = min(E) + Eint/3;
  HE = LE + Eint/3;
  
  %% evaluating and buying process
  for i = 1:simtime
    fund = fund + fb(i); % add this yr's fund to the account
    %[cumb,fund,buy] = buystrat(buy,code(8),cumb,fund,fb(i),ben(i),C(i),E(i),ff(i),fr(i),al,be,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf,LE,HE);
    [cumb,fund,buy] = buystrat(buy,code(mod(god,length(code))+1),cumb,fund,fb(i),ben(i),C(i),E(i),ff(i),fr(i),al,be,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf,LE,HE);
    fundt(i) = fund;
  end
  strcumb(mod(god,length(code))+1) = strcumb(mod(god,length(code))+1) + cumb;
  %fprintf('cumb=%.2f\n',cumb);
  fprintf('bought %d times\n',sum(buy));
  cumb = 0;
  %% calculate correlation btw f and C,ben, tjs.
  %fprintf("corelations between f and C,B,tj\n");
  cor = corrcoef(f(1:simtime),C);
  mfCcor = mfCcor + cor(1,2);
  %fprintf('cor(f,C)=%.2f\n',cor(1,2));
  cor = corrcoef(f(1:simtime),ben);
  mfBcor = mfBcor + cor(1,2);
  %fprintf('cor((f(1:simtime),tjs);
  cor = corrcoef(f(1:simtime),tjs);
  mftjcor = mftjcor + cor(1,2);
  %fprintf('cor(f,tj)=%.2f\n',cor(1,2));
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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PLOTTING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  plotting = 0;
  plw = [1,50];
  if plotting == 1
    what2pl = [0,11,9]; % [5,6,7];
    tiledlayout(length(what2pl),1)
    for pl = 1:length(what2pl)
      if what2pl(pl) == 0 % ff, fr
        nexttile
        plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
        xlabel('time');
        ylabel('net return');
        hold on
        plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
        legend('ff','fr');
        hold off
      elseif what2pl(pl) == 1 %f, ff, fr, fb
        nexttile
        plot(t(plw(1):plw(2)),f(plw(1):plw(2)),'k');
        xlabel('time');
        ylabel('net return');
        hold on
        plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
        plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
        plot(t(plw(1):plw(2)),fb(plw(1):plw(2)),'b');
        legend('f','ff','fr','fb');
        hold off
      elseif what2pl(pl) == 2 % fr-ff & f
        nexttile
        plot(t(plw(1):plw(2)),(fr(plw(1):plw(2))-ff(plw(1):plw(2))));
        %plot(t(plw(1):plw(2)),f(plw(1):plw(2)),'k');
        xlabel('time');
        ylabel('net return');

      elseif what2pl(pl) == 3 % f
        nexttile 
        plot(plw(1):plw(2),f(plw(1):plw(2)),'k');
        legend('f');
      elseif what2pl(pl) == 4 % c
        %land cost
        nexttile
        plot(t(plw(1):plw(2)),C);%.2f\n',cor(1,2));
        legend('land cost');
      elseif what2pl(pl) == 5 % ff
        nexttile
        plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
        legend('ff');
      elseif what2pl(pl) == 6 % fr
        nexttile
        plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
        legend('fr');
      elseif what2pl(pl) == 7 % fb
        nexttile
        plot(t(plw(1):plw(2)),fb(plw(1):plw(2)),'r');
        legend('fb');
      elseif what2pl(pl) == 8 % fund
        nexttile 
        plot(t(plw(1):plw(2)),fundt);
        legend('fund');
        xlabel('time');
      elseif what2pl(pl) == 9 % tj over time
        nexttile
        plot(t(plw(1):plw(2)),tjs(plw(1):plw(2)))
        legend('tj')
        xlabel('time')
      elseif what2pl(pl) == 10 % tj histogram
        nexttile
        histogram(tjs);
        xlabel('tj');
      elseif what2pl(pl) == 11 % C histogram
        nexttile
        histogram(C(C~=0)); % used when plotting only distribution with certain t_j
        Ctitle = strcat(Ctitle, sprintf('af=%.1f,ar=%.1f, tj%d fr',af,ar,cind));
        title(Ctitle)
        %histogram(C);
        xlabel('C');
      elseif what2pl(pl) == 12  % buying points
        nexttile
        scatter(t(plw(1):plw(2)),buy);
        legend('buying points')
      else  % B histogram
        nexttile
        histogram(ben);
        xlabel('ben');
      end
    end
  end
  if mod(god,100) == 0
    fprintf('god=%d\n',god);
  end
end

%tjsrecep = tjsrecep/(godsim/length(code));
%plot(1:length(tjsrecep),tjsrecep);
%xlabel('tjsrecep');
strcumb = strcumb/(godsim/length(code));
receptacle = {strcumb,t,ff,fr,fb,tjs,C,fundt,ben,buy};
%fprintf("str     mean cumb\n");
%stratstr = {'CVAL','Lc','Hc','Lff','Hff','Lfr','Hfr','LE','HE'};
%for i = 1:length(strcumb)
%  fprintf("%s       %.2f\n",stratstr{i},strcumb(i));
%end


% fprintf('mean of cor(f,C)=%.2f\n',mfCcor/godsim);
% fprintf('mean of cor(f,B)=%.2f\n',mfBcor/godsim);
% fprintf('mean of cor(f,tj)=%.2f\n',mftjcor/godsim);
%fprintf("mean tj=%.2f\n",sum(tjs)/simtime);
%fprintf("pr(tj=1) = %.2f\n", length(tjs(tjs==1))/length(tjs));
%fprintf("pr(tj=2) = %.2f\n", length(tjs(tjs==2))/length(tjs));
%fprintf("pr(tj=3) = %.2f\n", length(tjs(tjs==3))/length(tjs));