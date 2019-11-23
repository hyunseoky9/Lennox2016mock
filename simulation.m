% simulation script. Uses function opt_act
% output of the script is the mean conservation value for borrowing and fixed model
% set parameter values and RUN!

%%%%%%% parameters to set%%%%%%%%%%
b_default = 5; % budget
z = 1; % conservation value parameter (conservation value) = (conservation status)^z
c_max = 10; % maximum cost of a parcel allowed (in $100,000 in the paper)
s_max = 10; % maximum cons. status of a parcel allowed
d_max = b; % maximum debt allowed to accrue. set to 0 if no borrowing allowed (fixed)
p = 0.2; % probability parameter in the binomial distribution that chooses parcel cost
n = 10; % 240 in the paper
rep = 1000; %  simulation rep num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modval = [0 0]; % model value calculator at the end
for mo = 1:2 % run sim 2 times with dif model
	if mo == 1 % fixed
		d_max = 0;
	else % borrow
		d_max = b;
	end
	optm = opt_act(b_default,z,c_max,s_max,d_max,p,n); % get best action space
	modval_accum = 0;
	for i = 1:rep % simulation start
		modvalpt = 0; % model value calculator per time
		d = 0; % initial debt 0
		for t = 1:n
			c = binornd(c_max-1,p) + 1; % c instance
			g = c*(c_max-1)/c_max^2; % binom prob parameter for s (dep. on c)
			s = binornd(s_max-1,g) + 1; % s instance
			b = b_default; % budget
			repay = ceil(d/3); % payment of the debt in the beginning of time step
			d = d - repay; % debt reduced
			if optm(s,c,d+1,t) == 1 % buy
				modvalpt = modvalpt + s^z; % add conservation value by buying
				if c > (b - repay) % debt increased if borrowed
					d = c - (b - repay);
				end
			end
		end
		modval_accum = modval_accum + modvalpt;
	end
	modval(mo) = modval_accum/rep; % mean overall value throughout time recorded per model
end

% print out results
fprintf("\n\nmean conservation value by model\n");
fprintf("fixed = %.3f\n",modval(1));
fprintf("borrow = %.3f\n",modval(2));