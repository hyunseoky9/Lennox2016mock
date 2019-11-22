%%command script for running sdp function
%% create optm tensor that stores all the information for best action given state variable
%% at certain timestep (size = c_max*s_max*n)
b = 5; % budget
z = 1; % conservation value parameter (conservation value) = (conservation status)^z
c_max = 10; % maximum cost of a parcel allowed (in $100,000 in the paper)
s_max = 10; % maximum cons. status of a parcel allowed
d_max = b; % maximum debt allowed to accrue. set to 0 if no borrowing allowed (fixed)
p = 0.2; % probability parameter in the binomial distribution that chooses parcel cost
n = 10; % 240 in the paper
probm = zeros([s_max,c_max]); % pre-calculation of probabilities of parcel character probm(2,3) = Pr(s=2,c=3)
term = 0;
for i = 1:s_max % calculating probm matrix
    for j = 1:c_max
        g = j*(c_max-1)/c_max^2; % g is the probability parameter in the binom. dist. that picks s value.(dep. on c)
        term = nchoosek(s_max-1,i-1)*g^(i-1)*(1-g)^(s_max-i)*...
            nchoosek(c_max-1,j-1)*p^(j-1)*(1-p)^(c_max-j);
        probm(i,j) = term;
    end
end
optm = zeros([s_max,c_max,d_max+1,n]); %optimal action table given state vector and time
bestvalarray = zeros([s_max,c_max,d_max+1]); % best value obtainable at a given time step
bestvalarray_new = zeros([s_max,c_max,d_max+1]);
tic % timer
% calculation of optm through backwards recursion
for k = n:-1:1 
	for m = 0:d_max % element-wise multiplication between probm and bestval matrix with fixed debt column.
		bestvalarray(:,:,m+1) = probm.*bestvalarray(:,:,m+1);
	end
	%if k==8
	%	for aa = 0:d_max
	%		disp(sum(bestvalarray(:,:,aa+1),'all'));
	%	end
	%	%disp(bestvalarray(:,:,:))
	%end
	for i = 1:s_max
		for j = 1:c_max
			for l = 0:d_max
				repay = ceil(l/3); % repayment for the debt carried over from last step
				netbudget = b - repay; % net budget after paying repayment
				if netbudget + (d_max - (l-repay)) >= j % if you can buy with borrowing
					if netbudget >= j % calculation for probability different depending on borrowing or not
						% next step debt variable will only have left over debt
						ss = sum(bestvalarray(:,:,l-repay+1),'all'); 
					else
						% next step d will have left over + borrowed $ this step.
						ss = sum(bestvalarray(:,:,l-repay+(j-netbudget)+1),'all'); 
					end
					[bestvalarray_new(i,j,l+1),optm(i,j,l+1,k)] = findbest(i,[0,1],z,ss);
				else % if you can't buy even with borrowing maximum amount allowed
					ss = sum(bestvalarray(:,:,l-repay+1),'all');
					[bestvalarray_new(i,j,l+1),optm(i,j,l+1,k)] = findbest(i,[0],z,ss);
				end
			end
		end
	end
%	if k==9
%		disp(bestvalarray_new)
%	end
	bestvalarray = bestvalarray_new;
end
toc % timer end
%write out optm into csv file

%% can later develop more to output the optimal action into csv 
%% for now computation time is super short so no need.
%fileID = fopen('optimal_aciton.csv','w');
%for k = 1:n
%	for j=1:c_max
%
%		fprintf(fileID,'a');
%	end
%end
%fclose(fileID);