function interestp = interested_parameters(parambundle)
% finds parameters that are arrays (which corresponds to parameters that I want to keep track of in the file)
% and collect their index
interestp = [];
for i = 1:length(parambundle)
  if length(parambundle{i}) >= 2
    interestp = [interestp i];
  end
end 