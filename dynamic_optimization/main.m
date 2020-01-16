clear all
rng(1);
pw= [600,625]; %plot window
t = linspace(1,800,800);

%general economy
a = 0.8;
f0r = 30;
f0 = f0r*(1-a);
sig = 10*(1-a^2);

%forestry
lf = 0.7; %lambda
af = 0.4;
ff0r = 25;
ff0 = ff0r*(1-lf*af);
sigf = 4*(1-lf^2*af^2);
factorf = lf*(1-af)*f0r/(1-af*lf);
%housing
lr = 0.7;
ar = 0.1;
fr0r = 25;
fr0 = fr0r*(1-lr*ar);
sigr = 2*(1-lr^2*ar^2);

%donation
lb = 0.7;
ab = 0.9;
fb0r = 25;
fb0 = fb0r*(1-lb*ab);
sigb = 4*(1-lb^2*ab^2);

for j = 1
  f = zeros(0,length(t));
  for  i = 1:length(t)
    if i == 1
      f(i) = (f0 + normrnd(0,1));
    else
      f(i) = a*f(i-1) + (f0 + sig*normrnd(0,1));
    end
  end
  ff = zeros(1,length(t));
  fr = zeros(1,length(t));
  fb = zeros(1,length(t));
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
  plot(t(pw(1):pw(2)),f(pw(1):pw(2)));
  xlabel('time');
  ylabel('net return');
  hold on
  plot(t(pw(1):pw(2)),ff(pw(1):pw(2)));
  plot(t(pw(1):pw(2)),fr(pw(1):pw(2)));
  plot(t(pw(1):pw(2)),fb(pw(1):pw(2)));
  legend('f','ff','fr','fb');
  hold off
end

%fprintf('\n');
%fprintf('analytic ff mean = %f\n', ff0r) %analytic ff mean
%fprintf('simulation ff mean = %f\n',mean(ff(200:800))) %simulation mean

%fprintf('analytic f mean = %f\n',f0r) %analytic ff mean
%fprintf('simulation f mean = %f\n',mean(f(200:800))) %simulation mean
%fprintf('\n');
cor = corrcoef(ff,f);
%fprintf("cor(ff,f)= %.2f\n",cor(1,2))
cor = corrcoef(fr,f);
%fprintf("cor(fr,f)= %.2f\n",cor(1,2))
cor = corrcoef(fb,f);
%fprintf("cor(fb,f)= %.2f\n",cor(1,2))

%% simulation
burnin = 201; % burn in first few values of net return as they have not converged yet
f = f(burnin:end);
ff = ff(burnin:end);
fb = fb(burnin:end);
fr = fr(burnin:end);

simtime = 50; % number of buying opportunities
fund = 0; % money saved
cumcval = 0; % cummulative conservation value


bfn = 1; % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2

rho = 0.1; % economic discount rate
del = 0.01; % ecological discount rate


efmu = 0; % mean of ef
efsig2 = 1; % var of ef

edisc = repelem(1+rho,length(f)).^(0:(length(f)-1)); %economic discount rate.
ecodisc = repelem(1+del,length(f)).^(0:(length(f)-1)); %ecological discount rate
al = 1;
be = 1;
threshold =0;

for i = 1:simtime
  fprintf('time step=%d\n',i);
  
  nf = f(i:end);
  nff = ff(i:end);
  nfb = fb(i:end);
  nfr = fr(i:end);
  %fprintf('length of fb=%d\n',size(fb,2));
  %fprintf('this yrs donation = %.2f\n',fb(1));
  fund = fund + nfb(1); % add this yr's fund to the account
  
  % getting e_fj,e_rj, and b (if bfn=2)
  mu = [efmu 10];
  sigma = [efsig2 0.8; 0.8 1];
  R = mvnrnd(mu,sigma,1);
  e_fj = R(1); % indiv forest var
  e_rj = normrnd(0,0); % indiv dev var

  if bfn == 1
    e_fj = normrnd(efmu,efsig2);
    b = 10;
  elseif bfn >= 2
    b = R(2); % for benefit
  end
  if bfn == 3
    b = b^(1/2);
  end
  f_fj = nff + e_fj;
  f_fj(f_fj<0) = 0;
  f_rj = nfr + e_rj;
  f_rj(f_rj<0) = 0;


  % getting clearing time t_j
  t_j = 1;
  nedisc = edisc(1:end-i+1); % new economic discount rate for this sim.
  tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./nedisc(1:(end-t_j+1)));
  while tjoptim <= 0
    tjoptim = tjoptim - (f_rj(t_j)-f_fj(t_j))/nedisc(t_j);
    t_j = t_j + 1;
  end

  if t_j > 1
      fprintf('tj was bigger than 1 on step %d; t_j=%d\n',i,t_j);
  end
  
  % cost
  if t_j == 1
    c = sum(f_rj(t_j:end)./nedisc(t_j:end));
  else
    c = sum(f_fj(1:(t_j-1))./nedisc(1:(t_j-1))) + sum(f_rj(t_j:end)./nedisc(t_j:end));
  end
  %v = f_rj(t_j:end)./edisc(t_j:end);
  %disp(v(1:15))
  %benefit
  necodisc = ecodisc(1:end-i+1); %ecological discount rate
  B = sum(b./necodisc(t_j:end)); % conservation value (incl' discount and threat component)
  %v = b./ecodisc(t_j:end);
  %disp(v(1:15));
  % evaluating and buying process
  cval = (B/c)^al*fund^be; % conservation value
  if cval >= threshold && c <= fund % buy
    cumcval = cumcval + B;
    fund = fund - c;
    fprintf('bought!\n');
    fprintf('remaining fund=%.2f\n',fund);
  else
    fprintf('cost=%.2f',c);
    fprintf('remaining fund=%.2f\n',fund);
  end
end


fprintf('cumcval = %.2f\n',cumcval);
