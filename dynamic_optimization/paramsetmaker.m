function [paramset] = paramsetmaker(parambundle)
% makes set of parameter vector with all the element combinations.
paramset = {};
paramsetnum  = 1;

arraycellidx = [];
param = [];
for i = 1:length(parambundle)
	if length(parambundle{i}) > 1
		paramsetnum = paramsetnum*length(parambundle{i});
		param = [param 0];
		arraycellidx = [arraycellidx i];
	else
    val = parambundle{i};
		param(end+1) = val;
	end
end



for i = 1:paramsetnum
  paramset{i} = param;
end

if isempty(arraycellidx)
  paramset{1} = param;
else
  arrays = {};
  for i = 1:length(arraycellidx)
    arrays{i} = parambundle{arraycellidx(i)};
  end
  pset = recursion({},1,length(arraycellidx),arrays,[]);
  for i = 1:paramsetnum
    param(arraycellidx) = pset{i};
    paramset{i} = param;
  end
end



  function pset = recursion(pset,idx,finish,arrays,set)
    if idx == finish
      for i = 1:length(arrays)
        nset = [set arrays(i)];
        pset{end+1} = nset;
      end
    else
      for i = 1:length(arrays{idx})
        array = arrays{idx};
        nset = [set array(i)];
        pset = recursion(pset,idx+1,finish,arrays{2:end},nset);
      end
    end
  end
end