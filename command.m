%%command script for running sdp function
b = 5;
x = 0;
y = [0,1];
z = 2;
c_max = 10;
s_max = 10;
p = 0.2;
n = 1; % 20 in the paper
r = 0;
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
initbval = zeros([s_max,c_max]); %best value with given initial state variables element i,j is for s=i c=j
initbact = zeros([s_max,c_max]); %same thing except for best action.
tic
[bval,bact] = sdp(b,x,y,r,z,s_max,c_max,probm,p,n);
toc
fprintf('n=%d\n',n);
fprintf('bval=%.5f, bact=%d\n\n\n',bval,bact);