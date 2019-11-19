function [bestval,bestaction] = sdp(b,x,y,r,p,n) 
%% function for stochastic dynamic programming
% input:
% b = budget; scalar
% x = state space; 1x2 array 
% y = action space (0=no buy, 1= buy); 1x2 array
% r = sum of all the possible probabilities of buying land with s=i, c=j;
% p = parameter vector; [z,s_max,c_max,p]
% output:
% n = number of times game is played
% val = expected return by choosing an action at initial time; 1x2 array
bestaction = 0;
bestval = 0;
val = zeros([1,length(y)]);
if n > 1
    for action = 1:length(y)
        if y(action) == 0 %don't buy
            val(action) = sdp(b,x,y,r,n-1);
        elseif y(action) == 1 % buy
            val(action) = 0;
            for i = 1:p(2)
                for j = 1:p(3)
                    g = j*(p(3)-1)/p(3)^2; %probability for conservation status determined by cost
                    term = nchoosxe(p(2)-1,i-1)*g^(i-1)*(1-g)^(p(2)-i)*...
                        nchoose(p(3)-1,j-1)*p(4)^(j-1)*(1-p(4))^(p(3)-j)/r;
                    val(action) = val(action) + term*(i^p(1) + sdp(b,x+i^p(1),y,r,p,n-1));
                end
            end
        else % borrow (for later)
        end
    end
else 
    for action = 1:length(y)
        if y(action) == 0
            val(action) = x;
        elseif y(action) == 1
            val(action) = x + r;
        else %borrow (for later)
        end
    end
end
[bestval,bestaction] = max(val);
bestaction = y(bestaction);

