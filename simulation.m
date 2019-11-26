% simulation script. Uses function opt_act
% output of the script is the mean conservation value for borrowing and fixed model
% set parameter values and RUN!

%%%%%%% parameters to set%%%%%%%%%%
b_default = 5; % budget
z = 2; % conservation value parameter (conservation value) = (conservation status)^z
c_max = 10; % maximum cost of a parcel allowed (in $100,000 in the paper)
s_max = 10; % maximum cons. status of a parcel allowed
d_max = b_default; % maximum debt allowed to accrue. set to 0 if no borrowing allowed (fixed)
p = 0.8; % probability parameter in the binomial distribution that chooses parcel cost
n = 240; % 240 in the paper
rep = 1000; %  simulation rep num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modvalm = [0 0]; % model value mean calculator at the end
modvalstd = [0 0]; % model value std at the end
for mo = 2:2 % run sim 2 times with dif model
	if mo == 1 % fixed
		d_max = 0;
	else % borrow
		d_max = b;
	end
	optm = opt_act(b_default,z,c_max,s_max,d_max,p,n); % get best action space
	modval_accum = zeros([1000,1]);
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

