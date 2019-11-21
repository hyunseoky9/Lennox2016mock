%%WRONG
function [bestval,bestaction] = sdp(b,x,y,r,z,s_max,c_max,probm,p,n) 
%% function for stochastic dynamic programming
% input:
% b = budget; scalar
% x = state space; 1x2 array 
% y = action space (0=no buy, 1= buy); 1x2 array
% r = sum of all the possible probabilities of buying land with s=i, c=j;
% s_max = maximum conservation status value
% c_max = maximum land cost value
% probm = probability matrix where probm(i,j) is Pr of s=i,c=j
% p = probability for binom dist for cost
% output:
% n = number of times game is played
% val = expected return by choosing an action at initial time; 1x2 array
bestaction = 0;
bestval = 0;
val = zeros([1,length(y)]);
if n > 1
    for action = 1:length(y)
        if y(action) == 0 %don't buy
            val(action) = sdp(b,x,y,r,z,s_max,c_max,probm,p,n-1); %%%%%
        elseif y(action) == 1 % buy
            val(action) = 0;
            for i = 1:s_max
                for j = 1:b
                    g = j*(c_max-1)/c_max^2; %probability for conservation status determined by cost
                    val(action) = val(action) + probm(i,j)*(i^z + sdp(b,x+i^z,y,r,z,s_max,c_max,probm,p,n-1));
                end
            end
        else % borrow (for later)
        end
    end
else %last recursion
    for action = 1:length(y)
        if y(action) == 0
            val(action) = x;
        elseif y(action) == 1
            for i = 1:s_max
                for j = 1:b
                    g = j*(c_max-1)/c_max^2; %probability for conservation status determined by cost
                    val(action) = val(action) + probm(i,j)*(i^z);
                end
            end
            val(action) = val(action) + x;
        else %borrow (for later)
        end
    end
end
[bestval,bestaction] = max(val);
bestaction = y(bestaction);







