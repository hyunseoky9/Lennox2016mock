function duh = pl(what2pl,tjv,n,receptacle,param,i)
% laying out plots

if any(what2pl == 12)
  m = length(what2pl) - 1 + length(tjv);
else
  m = length(what2pl);
end

for j = 1:length(what2pl)
  if what2pl(j) == 12 % if contingent c dist has to be plotted for multiple contingents
    for k = 1:length(tjv)
      subplot(m,n,(n*k-(n-1))+(i-1))
      bleh = plotting(receptacle,param,what2pl(j),tjv(k));
    end
  else
    subplot(m,n,(n*j-(n-1))+(i-1))
    bleh = plotting(receptacle,param,what2pl(j),0); %make it better
  end
end
