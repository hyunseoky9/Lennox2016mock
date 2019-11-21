%%command script for running sdp function
b = 5;
z = 1;
c_max = 10;
s_max = 10;
p = 0.2;
n = 1; % 240 in the paper
probm = zeros([s_max,c_max]);
term = 0;
for i = 1:s_max
    for j = 1:c_max
        g = j*(c_max-1)/c_max^2;
        term = nchoosek(s_max-1,i-1)*g^(i-1)*(1-g)^(s_max-i)*...
            nchoosek(c_max-1,j-1)*p^(j-1)*(1-p)^(c_max-j);
        probm(i,j) = term;
    end
end
optm = zeros([s_max,c_max,n]); %optimal action table given state vector and time
bestvalarray = zeros([s_max,c_max]);
bestvalarray_new = zeros([s_max,c_max]);
tic
for k = n:1
	ss = sum(probm.*bestvalarray,'all');
	for i = 1:s_max
		for j = 1:c_max
			if j > b
				[bestvalarray_new(i,j),optm(i,j,k)] = findbest(i,j,d,[0],z,ss,k);
			else
				[bestvalarray_new(i,j),optm(i,j,k)] = findbest(i,j,d,[0,1],z,ss,k);
			end
		end
	end
	bestvalarray = bestvalarray_new;
end
toc