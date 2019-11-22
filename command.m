%%command script for running sdp function
%% create optm tensor that stores all the information for best action given state variable
%% at certain timestep (size = c_max*s_max*n)
b = 5;
z = 1;
c_max = 10;
s_max = 10;
d_max = b;% set to 0 if no borrowing allowed (fixed)
p = 0.2;
n = 240; % 240 in the paper

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
optm = zeros([s_max,c_max,d_max+1,n]); %optimal action table given state vector and time
bestvalarray = zeros([s_max,c_max,d_max+1]);
bestvalarray_new = zeros([s_max,c_max,d_max+1]);
tic
for k = n:-1:1
	ss = sum(probm.*bestvalarray(:,:,SOMETHING_HERE),'all');
	for i = 1:s_max
		for j = 1:c_max
			for l = 0:d_max
				repay = ceil(l/3);
				netbudget = b - repay;
				if netbudget + d_max >= c % can buy
					[bestvalarray_new(i,j,l-repay),optm(i,j,l-repay,k)] = findbest(i,j,l-repay,[0,1],z,ss,k,);
				elseif netbudget < c
					[bestvalarray_new(i,j,l-repay),optm(i,j,l-repay,k)] = findbest(i,j,l,[0],z,ss,k);
				end
			end
		end
	end
	bestvalarray = bestvalarray_new;
end
toc
%disp(bestvalarray);
%disp(optm);
%write out optm into csv file
fileID = fopen('optimal_aciton.csv','w');
for k = 1:n
	for j=1:c_max

		fprintf(fileID,'a');
	end
end

fclose(fileID);