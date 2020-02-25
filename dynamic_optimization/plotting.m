function [foo] = plotting(receptacle,param,what2pl,tjv)
strcumb = receptacle{1};
t = receptacle{2};
x = receptacle{3};
xf = receptacle{4};
xr = receptacle{5};
xb = receptacle{6};
tjs = receptacle{7};
C = receptacle{8};
fundt = receptacle{9};
ben = receptacle{10};
buy = receptacle{11};

%tjv = [1,5,9,13,15];
simtime = param(34);
plw = [100,400];
%what2pl = [12]; % [5,6,7];

%% making tile layout.
%if length(what2pl) <= 5 && any(what2pl == 12) ~= 1
%  tiledlayout(length(what2pl),1)
%else
%  if any(what2pl == 12)
%    pln = length(what2pl) - 1 + length(tjv); % plot number
%  else 
%    pln = length(what2pl);
%  end
%  tiledlayout(5,ceil(pln/5))
%end

for pl = 1:length(what2pl)
  if what2pl(pl) == 0 % xf, xr
    %nexttile
    plot(t(plw(1):plw(2)),xf(plw(1):plw(2)),'g');
    xlabel('time');
    ylabel('net return');
    hold on
    plot(t(plw(1):plw(2)),xr(plw(1):plw(2)),'r');
    legend('xf','xr');
    hold off
  elseif what2pl(pl) == 1 %x, xf, xr, xb
    %nexttile
    plot(t(plw(1):plw(2)),x(plw(1):plw(2)),'k');
    xlabel('time');
    ylabel('net return');
    hold on
    plot(t(plw(1):plw(2)),xf(plw(1):plw(2)),'g');
    plot(t(plw(1):plw(2)),xr(plw(1):plw(2)),'r');
    plot(t(plw(1):plw(2)),xb(plw(1):plw(2)),'b');
    legend('x','xf','xr','xb');
    hold off
  elseif what2pl(pl) == 2 % xr-xf & x
    %nexttile
    plot(t(plw(1):plw(2)),(xr(plw(1):plw(2))-xf(plw(1):plw(2))));
    %plot(t(plw(1):plw(2)),x(plw(1):plw(2)),'k');
    xlabel('time');
    ylabel('net return');

  elseif what2pl(pl) == 3 % x
    %nexttile 
    plot(plw(1):plw(2),x(plw(1):plw(2)),'k');
    legend('f');
  elseif what2pl(pl) == 4 % c
    %land cost
    %nexttile
    plot(t(plw(1):plw(2)),C(plw(1):plw(2)));%.2f\n',cor(1,2));
    legend('land cost');
  elseif what2pl(pl) == 5 % xf
    %nexttile
    plot(t(plw(1):plw(2)),xf(plw(1):plw(2)),'g');
    legend('xf');
  elseif what2pl(pl) == 6 % xr
    %nexttile
    plot(t(plw(1):plw(2)),xr(plw(1):plw(2)),'r');
    legend('xr');
  elseif what2pl(pl) == 7 % xb
    %nexttile
    plot(t(plw(1):plw(2)),xb(plw(1):plw(2)),'r');
    legend('xb');
  elseif what2pl(pl) == 8 % fund
    %nexttile 
    plot(t(plw(1):plw(2)),fundt(plw(1):plw(2)));
    legend('fund');
    xlabel('time');
    ylim([0 5000]);
  elseif what2pl(pl) == 9 % tj over time
    %nexttile
    plot(t(plw(1):plw(2)),tjs(plw(1):plw(2)))
    legend('tj')
    xlabel('time')
  elseif what2pl(pl) == 10 % tj histogram
    %nexttile
    histogram(tjs);
    xlabel('tj');
    title('tj histogram')
  elseif what2pl(pl) == 11 % C histogram
    %nexttile
    histogram(C); % used when plotting only distribution with certain t_j
    fprintf('length(C)=%d\n',length(C));
    xlim([200,400]);
    ylim([0 1000]);
    %histogram(C);
    xlabel('C');
    title("C histogram");
  elseif what2pl(pl) == 12 % C histogram specific contingent
    %nexttile
    val = tjv;
    I = find(tjs == val);
    histogram(C); % used when plotting only distribution with certain t_j
    histogram(C(I)); % used when plotting only distribution with certain t_j
    Ctitle = '';
    Ctitle = strcat(Ctitle, sprintf('af=%.1f,ar=%.1f, tj%d',param(6),param(10),val));
    title(Ctitle)
    xlim([240,350]);
    %ylim([0 620]);
    %histogram(C);
    xlabel('C');
  elseif what2pl(pl) == 13  % buying points
    %nexttile
    scatter(t(plw(1):plw(2)),buy(plw(1):plw(2)));
    legend('buying points')
  elseif what2pl(pl) == 14 % B histogram
    %nexttile
    histogram(ben);
    xlabel('ben');
    title('histogram of B')
    [N,edges] = histcounts(ben);
  elseif what2pl(pl) == 15 % histogram of ecological benefits of bought lands
    %nexttile
    I = find(buy==1);
    histogram(ben(I));
    xlabel('ben')
    title('histogram of eco benefit of bought lands');
  elseif what2pl(pl) == 16 % mean of xr's 10 steps ahead
  	%nexttile
  	xrm = zeros(1,(plw(2)-plw(1)+1));
  	for i = plw(1):plw(2)
  		xrm(i) = mean(xr(i:(i+10)));
  	end
  	plot(t(plw(1):plw(2)),xrm(plw(1):plw(2)));
  	legend('xr fw mean');
    ylim([22 32]);
  	xlabel('time')
  elseif what2pl(pl) == 17
    plot(C,ben);
    words = sprintf('cost and benefit plot corr=%.2f',corrcoef(C,ben));
    title(words);
  else % plot of accumulation of B of bought lands in descending order 
    %nexttile
    I = find(buy==1);
    sortben = sort(ben(I),'descend');
    cum = cumsum(sortben);
    upto = 100;
    plot((1:upto),cum(1:upto));
    title('accumulation of B of bought lands in descending order')
  end
end

foo = 1;