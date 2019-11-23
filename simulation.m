b = 5; % budget
z = 1; % conservation value parameter (conservation value) = (conservation status)^z
c_max = 10; % maximum cost of a parcel allowed (in $100,000 in the paper)
s_max = 10; % maximum cons. status of a parcel allowed
d_max = b; % maximum debt allowed to accrue. set to 0 if no borrowing allowed (fixed)
p = 0.2; % probability parameter in the binomial distribution that chooses parcel cost
n = 10; % 240 in the paper
for mo = 1:2 % run sim 2 times with dif model
	if mo == 1 % fixed
		d_max = 0;
	else % borrow
		d_max = b;
	end
end

optm = opt_act(b,z,c_max,s_max,d_max,p,n);

rep = 1000; % repetition simulation
val_f = 0; % consrevation value counter for fixed
val_b = 0; % counter for borrow
for i = 1:rep
	for t = 1:n
		c = binornd(c_max,p);
		g = c*(c_max-1)/c_max^2;
		s = binornd(s_max,g);
		val_f = optm()
	end
end
