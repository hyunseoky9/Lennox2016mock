function [bestval,bestaction] = findbest(b,c,s,d,y,z,s_max,c_max,probm,n,barray)
%% update optimal action matrix and best value for the given state and time.
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
% barrray = matrix of best values at time n+1 i,j is best value at s=i c=j
% output:
% bestval = best value given best action
% bestaction = best action given state 
val = zeros([1,length(y)])
bestaction = 0;
bestval = 0;
if n ~= 1 %  recursion
	for action = 1:length(y)
		if y(action) == 0 % don't buy/can't buy
			val(action) = 0 + ;
		else if y(action) == 1 % buy
			val(action) = s^z + probm.*barray;
		else % borrow and buy
			foo = 1;
			%do it later

else %last step in recursion
	for action = 1:length(y)
		if y(action) == 0
			val(action) = 0
end
[bestval,bestaction] = max(action);
bestaction = y(bestaction);

