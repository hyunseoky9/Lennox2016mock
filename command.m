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
for k = n:1
	for i = 1:s_max
		for j = 1:c_max
			if j > b
				[bestvalarray_new(i,j),optm(i,j,k)] = findbest(b,j,i,d,[0],z,s_max,c_max,probm,n,bestvalarray);
			else
				[bestvalarray_new(i,j),optm(i,j,k)] = findbest(b,j,i,d,[0,1],z,s_max,c_max,probm,n,bestvalarray);
			end
		end
	end
end


tic
[bval,bact] = sdp(b,x,y,r,z,s_max,c_max,probm,p,n);
toc
fprintf('n=%d\n',n);
fprintf('bval=%.5f, bact=%d\n\n\n',bval,bact);