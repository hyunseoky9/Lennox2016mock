clear all
rng(1);
pw= [600,625]; %plot window
<<<<<<< HEAD
t = linspace(1,800,800); 
=======
t = linspace(1,800,800); % time 


>>>>>>> test
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

fprintf('\n');
fprintf('analytic ff mean = %f\n', ff0r) %analytic ff mean
fprintf('simulation ff mean = %f\n',mean(ff(200:800))) %simulation mean

fprintf('analytic f mean = %f\n',f0r) %analytic ff mean
fprintf('simulation f mean = %f\n',mean(f(200:800))) %simulation mean
fprintf('\n');
cor = corrcoef(ff,f);
fprintf("cor(ff,f)= %.2f\n",cor(1,2))
cor = corrcoef(fr,f);
fprintf("cor(fr,f)= %.2f\n",cor(1,2))
cor = corrcoef(fb,f);
fprintf("cor(fb,f)= %.2f\n",cor(1,2))

%% simulation
burnin = 200; % burn in first few values of net return as they have not converged yet
f = f(burnin:end);
ff = ff(burnin:end);
fb = fb(burnin:end);
fr = fr(burnin:end);

simtime = 1; % number of buying opportunities
fund = 0; % money saved
cumcval = 0;
for i = 1:simtime
  f = f(i:end);
  ff = ff(i:end);
  fb = fb(i:end);
  fr = fr(i:end);

  bfn = 1; % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2

  e_rj = normrnd(0,0); % indiv dev var
  efmu = 0;
  efsig = 1;
  if bfn == 1
    e_fj = normrnd(efmu,efsig);
    b = 10;
  else
    mu = [efmu 10];
    sigma = [efsig 1.5; 1.5 1];
    R = mvnrnd(mu,sigma,1);
    e_fj = R(1); % indiv forest var
    b = R(2); % for benefit
  end
  f_fj = ff _ e_fj;
  f_fj(f_fj<0) = 0;
  f_rj = fr _ e_rj;
  f_rj(f_rj<0) = 0;

  rho = 0.1; % economic discount rate
  del = 0.01; % ecological discount rate

  % getting clearing time t_j
  t_j = 1;
  edisc = repelem(1+rho,length(f)-burnin+1).^(1:length(f)-burnin+1); %economic discount rate.
  tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./edisc(1:(end-t_j+1)));
  while tjval <= 0
    t_j = t_j + 1;
    tjoptim = sum((f_rj(t_j:end)-f_fj(t_j:end))./edisc(1:(end-t_j+1)));
  end

  % cost
  if t_j == 1
    c = sum(f_rj(t_j:end)./edisc(t_j:end));
  else
    c = sum(f_fj(1:(t_j-1))./edisc(1:(t_j-1))) + sum(f_rj(t_j:end)./edisc(t_j:end));
  end

  %benefit
  if bfn == 3
    b = b^(1/2);
  end
  ecodisc = repelem(1+del,length(f)-burnin+1).^(1:length(f)-burnin+1); %ecological discount rate
  B = sum(b./ecodisc(t_j:end)); % conservation value (incl' discount and threat component)

  % evaluating and buying process
  al = 1;
  be = 1;
  threshold =0;
  cval = (B/c)^al*fund^be; % conservation value
  if cval >= threshold % buy
    cumcval = cumcval + B;
  else % save
    fund = fund + fb(1);
  end
end
commit1 = 0;
commit2 = 1;
commit3 = 2;

whatsup = 0;
whatsup2 = 1;


fprintf('cumcval = %.2f',cumcval)';

