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
