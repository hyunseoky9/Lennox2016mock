
function [x,y] = autoregression(t,s,o,c,x0,y0)
%% simulation of real estate and forestry market value.
% input:
% t = number of timestep
% s = variance in epsilon factor
% o = order of AR
% c = coefficients for 
% x0 = initial x state for first timestep estimation
% y0 = same for y.
% output:
% x = real estate market value for all timestep
% y = forestry market value for all timestep

x = zeros([1 t+o]); % real estate market var
y = x; % timber market var
x(1:o) = x0;
y(1:o) = y0;
for i = 1:t
	x(i+o) = sum(c.*x((i-1+o):-1:i)) + normrnd(0,s); % i-o+o=i;
	y(i+o) = sum(c.*y((i-1+o):-1:i)) + normrnd(0,s);
end
%plot(x[(o+1):length(x)]~t,type='l',ylim=c(min(min(y),min(x)),max(max(y),max(x))))
%lines(y[(o+1):length(y)]~t,col='red')
  
