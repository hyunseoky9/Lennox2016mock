% t_j when overall real estate value is higher than overall forestry value. 
% different method of generating f_f f_r and b. 
% 1) current random autoregression one for forestry and another for real estate.
% 2) one autoregression for general economy and forestry, real estate, and budget are causal to the gen. economy.
% 3) make periodic function for forestry, real estate, and budget, and have a parameter that mess with their phase (review ode)
%		Also, make parameters for amplitude and cycle center for each value function.
%  1 not advisable. Let's go with 2 and 3.

%%%%%%%% parameters %%%%%%%%%
%autoregression params%%%
T = 200; % timestep
o = 3; % order
s = 0.5; %variance
% x0 and y0 have same mean
x0 = [12,14,13];
y0 = [13,14,12]; % initial 
coef = [0.33,0.33,0.34]; %coefficients for autoregression. sum is around 1 to keep the x or y from increasing or decreasing too much
%%%%%%%%%%%%%%%%%%%%%%%%%%
rep = 1000; %  simulation rep num
c = 0; % cost of a parcel
s = 0; % conservation status of a parcel
z = 1; % power factor for conservation status
d_max = 0; % maximum debt allowed to accrue. set to 0 if no borrowing allowed (fixed)
efu = 0; % mean of the error term for forest parcel value from the mean.
efsd = 0; % sd of the ""
eru = 0; % mean of the error term for real estate market parcel 
ersd = 0; % sd of the ""
be = 0; % budget error term. Budget is the error term + r.e. market price.
tj = 0; % conversion time of the parcel that's out for sell
numo = 5; % number of heuristics method
rho = 0; % economic discount rate 
del = 0; % ecological discount rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("AR params:\n");
fprintf("t=%d, o=%d, s=%.2f\n",T,o,s);
fprintf("AR coeffs = ")
disp(c);
fprintf("simulation parameters:\n");
fprintf("...");
modvalm = zeros([1 numo]); % model value mean calculator at the end
modvalstd = zeros([1 numo]; % model value std at the end
[x,y] = autoregression(t,s,o,coef,x0,y0);
for mo = 1:numo % run sim 2 times with dif model
	modval_accum = zeros([1000,1]);
	for i = 1:rep % simulation start
		modvalpt = 0; % model value calculator per time
		d = 0; % initial debt 0
		for t = 1:n
			ef = normrnd(efu,efsd);
			er = normrnd(eru,ersd);
			tj = ;% conversion time
			c = ;% something
			s = ef + 2*efsd;
			if s < 0
				s = 0;
			end
			be = normrnd(0,1);
			b = x + be; % budget
			repay = ceil(d/3); % payment of the debt in the beginning of time step
			d = d - repay; % debt reduced
			% have to specify if you can buy the parcel given net budget and debt you have (d<=d_max)
			if strategy(mo) == 1 % buy
				consv = ; % conservation value. function of s^z
				modvalpt = modvalpt + s^z; % add conservation value by buying
				if c > (b - repay) % debt increased if borrowed
					d = c - (b - repay);
				end
			end
		end
		modval_accum(i) = modvalpt;
	end
	modvalm(mo) = mean(modval_accum); % mean overall value throughout time recorded per model
	modvalstd(mo) = std(modval_accum);
end



% print out results
fprintf("\n\nmean and 1std confidence interval of conservation value by model\n");
fprintf("fixed mean= %.3f\n",modvalm(1));
fprintf("borrow mean= %.3f\n",modvalm(2));
fprintf("fixed C.I. = (%.3f, %.3f)\n",modvalm(1)-modvalstd(1),modvalm(1)+modvalstd(1));
fprintf("borrow C.I.= (%.3f, %.3f)\n",modvalm(2)-modvalstd(2),modvalm(2)+modvalstd(2));
