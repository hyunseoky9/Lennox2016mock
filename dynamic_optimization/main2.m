clear all
mftjcor = 0;
mfCcor = 0;
mfBcor = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rng(1);
timeline = 10000;
t = linspace(1,timeline,timeline);

%general economy
a = 0.8;
f0r = 30;
f0 = f0r*(1-a);
fsig2 = 10;
sig = fsig2*(1-a^2);

%forestry
lf = 0.7; %lambda
af = 0.4;
ff0r = 25;
ff0 = ff0r*(1-lf*af);
ffsig2 = 4;
sigf = ffsig2*(1-lf^2*af^2);
factorf = lf*(1-af)*f0r/(1-af*lf);

%housing
lr = 0.7;
ar = 0.4;
fr0r = 25;
fr0 = fr0r*(1-lr*ar);
frsig2 = 4;
sigr = frsig2*(1-lr^2*ar^2);

%donation
lb = 0.7;
ab = 0.9;
fb0r = 25;
fb0 = fb0r*(1-lb*ab);
fbsig2 = 4;
sigb = fbsig2*(1-lb^2*ab^2);

% buy strategy stuff
code = 1:9;
al = 1;
be = 1;

godsim = 1;%length(code)*2;
strcumb = zeros(1,length(code)); % strategy's cumulative benefit
for god = 1:godsim

  f = zeros(0,length(t));
  ff = zeros(1,length(t));
  fr = zeros(1,length(t));
  fb = zeros(1,length(t));
  regsin = 0; % 0=autoregression 1=sinusoidal (for ff,fr,fb)
  if regsin == 0   
    for  i = 1:length(t)
      if i == 1
        f(i) = (f0 + normrnd(0,1));
      else
        f(i) = a*f(i-1) + (f0 + sig*normrnd(0,1));
      end
    end

    for i = 1:length(t)
      if i == 1
        ff(i) = (ff0 + sigf*normrnd(0,1));
        fb(i) = (fb0 + sigb*normrnd(0,1));
        f(i) = (fr0 + sigr*normrnd(0,1));
      else
        ff(i) = lf*(af*ff(i-1) + (1-af)*(f(i) - f0r)) + (ff0 + sigf*normrnd(0,1));
        fb(i) = lb*(ab*fb(i-1) + (1-ab)*(f(i) - f0r)) + (fb0 + sigb*normrnd(0,1));
        fr(i) = lr*(ar*fr(i-1) + (1-ar)*(f(i) - f0r)) + (fr0 + sigr*normrnd(0,1));
      end
    end
  else
    period = 10;
    lag = 0;
    A = 1; % amplitude
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
  cor = corrcoef(ff,f);
  %fprintf("cor(ff,f)= %.2f\n",cor(1,2))
  cor = corrcoef(fr,f);
  %fprintf("cor(fr,f)= %.2f\n",cor(1,2))
  cor = corrcoef(fb,f);
  %fprintf("cor(fb,f)= %.2f\n",cor(1,2))
  %fprintf("\n");

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % simulation of tj, b, c, buying, etc.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  burnin = 201; % burn in first few values of net return as they have not converged yet
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
  cvalth = 80; % buying threshold for buystrat code 1

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

  edisc = repelem(1+rho,length(f)).^(0:(length(f)-1)); %economic discount rate.
  ecodisc = repelem(1+del,length(f)).^(0:(length(f)-1)); %ecological discount rate
  al = 1;
  be = 1;

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
   

    e_rj = normrnd(0,0); % indiv dev var

    if bfn == 1
      e_fj = normrnd(efmu,efsig2);
      b = 10;
    elseif bfn >= 2
      b_mu = 10;
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
    t_j = 1;
    nedisc = edisc(1:end-i+1); % new economic discount rate for this sim.
    tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./nedisc(1:(end-t_j+1))) - ch;
    while tjoptim <= 0 && t_j <= length(f_fj)
      tjoptim = tjoptim - (f_rj(t_j)-f_fj(t_j))/nedisc(t_j);
      t_j = t_j + 1;
    end

    tjs(i) = t_j;
    if t_j > 1
        %fprintf('tj was bigger than 1 on step %d; t_j=%d\n',i,t_j);
    end
    
    % cost
    if t_j == 1
      c = sum(f_rj(t_j:end)./nedisc(t_j:end));
    else
      c = sum(f_fj(1:(t_j-1))./nedisc(1:(t_j-1))) + sum(f_rj(t_j:end)./nedisc(t_j:end));
    end
    C(i) = c;
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
  cint = max(C)-min(C);
  Lc = min(C) + cint/3; % low threshold for cost
  Hc = Lc + cint/3; % high threshold for cost

  E = C./ben; % roi ratio
  Eint = max(E) - min(E);
  LE = min(E) + Eint/3;
  HE = LE + Eint/3;
  
  %% evaluating and buying process
  for i = 1:simtime
    fund = fund + fb(i); % add this yr's fund to the account
    [cumb,fund,buy] = buystrat(buy,code(8),cumb,fund,fb(i),ben(i),C(i),E(i),ff(i),fr(i),al,be,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf,LE,HE);
    %[cumb,fund,buy] = buystrat(buy,code(mod(god,length(code))+1),cumb,fund,fb(i),ben(i),C(i),E(i),ff(i),fr(i),al,be,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf,LE,HE);
  end
  strcumb(mod(god,length(code))+1) = strcumb(mod(god,length(code))+1) + cumb;
  fprintf("cumb=%.2f\n",cumb);
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

  %fprintf('cumb = %.2f\n',cumb);
  %fprintf('mincval = %.2f, maxcval = %.2f\n',mcval(1),mcval(2));
  %fprintf("mean cval = %.2f\n",cvalm/simtime);
  %fprintf('buy=%d \n',find(buy > 0));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PLOTTING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  plotting = 1;
  if plotting == 1
    what2pl = [0,5];
    tiledlayout(length(what2pl),1)
    for pl = 1:length(what2pl)
      if what2pl(pl) == 0 % ff, fr
        nexttile
        plot(t(1:simtime),ff(1:simtime));
        xlabel('time');
        ylabel('net return');
        hold on
        plot(t(1:simtime),fr(1:simtime));
        legend('ff','fr');
        hold off
      elseif what2pl(pl) == 1 %f, ff, fr, fb
        nexttile
        plot(t(1:simtime),f(1:simtime));
        xlabel('time');
        ylabel('net return');
        hold on
        plot(t(1:simtime),ff(1:simtime));
        plot(t(1:simtime),fr(1:simtime));
        plot(t(1:simtime),fb(1:simtime));
        legend('f','ff','fr','fb');
        hold off
      elseif what2pl(pl) == 2 % fr-ff & f
        nexttile
        plot(t(1:simtime),(fr(1:simtime)-ff(1:simtime)));
        plot(t(1:simtime),f(1:simtime));
        xlabel('time');
        ylabel('net return');

      elseif what2pl(pl) == 3 % f
        nexttile 
        plot(1:simtime,f(1:simtime));
        legend('f')
      elseif what2pl(pl) == 4 % c
        %land cost
        nexttile
        plot(t(1:simtime),C);%.2f\n',cor(1,2));
        legend('land cost')
      else 
        nexttile
        histogram(tjs);
        xlabel('tj');
      end
    end

    % buying points
    %nexttile
    %scatter(t(1:simtime),buy);
    %legend('buying points')

  end
  fprintf("god=%d\n",god);
end

strcumb = strcumb/(godsim/length(code));
fprintf("str     mean cumb\n");
stratstr = {'CVAL','Lc','Hc','Lff','Hff','Lfr','Hfr','LE','HE'};
for i = 1:length(strcumb)
  fprintf("%s       %.2f\n",stratstr{i},strcumb(i));
end
% fprintf('mean of cor(f,C)=%.2f\n',mfCcor/godsim);
% fprintf('mean of cor(f,B)=%.2f\n',mfBcor/godsim);
% fprintf('mean of cor(f,tj)=%.2f\n',mftjcor/godsim);
fprintf("mean tj=%.2f\n",sum(tjs)/simtime);
fprintf("pr(tj=1) = %.2f\n", length(tjs(tjs==1))/length(tjs));
fprintf("pr(tj=2) = %.2f\n", length(tjs(tjs==2))/length(tjs));
fprintf("pr(tj=3) = %.2f\n", length(tjs(tjs==3))/length(tjs));