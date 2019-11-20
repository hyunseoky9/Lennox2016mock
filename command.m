%%command script for running sdp function
b = 5;
x = 0;
y = [0,1];
z = 1;
c_max = 10;
s_max = 10;
p = 0.2;
n = 5; % 20 in the paper
r = 0;
for i = 1:10
    for j = 1:b
        g = j*(c_max-1)/c_max^2;
        r = r + nchoosek(s_max,i)*g^i*(1-g)^(s_max-i)*...
            nchoosek(c_max,j)*p^j*(1-p)^(c_max-j);
    end
end
[bval,bact] = sdp(b,x,y,r,z,s_max,c_max,p,n);
