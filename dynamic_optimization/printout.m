function printout(pname,param,interestp,stratstr,receptacle)
for i = interestp
  fprintf('%s=%.2f, ',pname{i},param(i));
end
fprintf("\n");
fprintf('str = %s\n',stratstr{param(end)})
fprintf('mean cumb:%.2f\n',receptacle{1});
fprintf('\n\n');
