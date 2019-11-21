function [bestval,bestaction] = brecursion(b,c,s,d,y,z,s_max,c_max,probm,n,fbestval)
%% bakcwards recursion method
% input:
% b = budget; scalar (state variable1)
% c = cost; scalar (state variable 2)
% s = conservation status; scalar (state variable 3)
% d = debt; scalar (state variable 4)
% y = action space; array (either {0},{0,1}, or {0,1,2}if there's borrowing)
% s_max = maximum conservation status value
% z = exponent for calculating conservation value from status.
% s_max = maximum status value
% c_max = maximum land cost value
% probm = probability matrix where probm(i,j) is Pr of s=i,c=j
% n = number of backwards recursions to do.
% fbestval = best value from forward time step
% p = probability for binom dist for cost
% output:
% bestval = maximum value for given set of choices
% bestaction = action choice in y that would yield maximum value.
val = zeros([1,length(y)])
bestaction = 0;
bestval = 0;
if n ~= 1 %  recursion
	for action = 1:length(y)
		if y(action) == 0 % don't buy/can't buy
			val(action) = 0;
		else if y(action) == 1 % buy
			val(action) = s^z;
		else % borrow and buy
			foo = 1;
			%do it later
else %last step in recursion
end

brecursion(b,c,s,d,y,z,s_max,c_max,probm,n-1,fbestval);
