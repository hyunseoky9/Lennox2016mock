function [bestval,bestaction] = findbest(c,s,d,y,z,ss,n)
%% update optimal action matrix and best value for the given state and time.
% input:
% c = cost; scalar (state variable 2)
% s = conservation status; scalar (state variable 3)
% d = debt; scalar (state variable 4)
% y = action space; array (either {0},{0,1}, or {0,1,2}if there's borrowing)
% s_max = maximum conservation status value
% z = exponent for calculating conservation value from status.
% ss = sum of element multiplication matrix of probm and besetvalarray matrix.
% n = number of backwards recursions to do.
% output:
% bestval = best value given best action
% bestaction = best action given state 
val = zeros([1,length(y)])
bestaction = 0;
bestval = 0;
for action = 1:length(y)
	if y(action) == 0 % don't buy/can't buy
		val(action) = 0 + ;
	else if y(action) == 1 % buy
		val(action) = s^z + ss;
	else % borrow and buy
		foo = 1;
		%do it later
end
[bestval,bestaction] = max(action);
bestaction = y(bestaction);
