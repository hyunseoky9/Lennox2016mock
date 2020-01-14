efmu = 0;
efsig = 1;
mu = [efmu 10];
sigma = [1 0.8; 0.8 1];
R = mvnrnd(mu,sigma,1);
e_fj = R(1); % indiv forest var
b = R(2); % for benefit
R