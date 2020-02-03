function [foo] = plotting(receptacle,param)
strcumb = receptacle{1};
t = receptacle{2};
ff = receptacle{3};
fr = receptacle{4};
fb = receptacle{5};
tjs = receptacle{6};
C = receptacle{7};
fundt = receptacle{8};
ben = receptacle{9};
buy = receptacle{10};
cind = receptacle{11};

simtime = param(34);
plw = [1,1000];
what2pl = [0,11,9,14]; % [5,6,7];
tiledlayout(length(what2pl),1)
for pl = 1:length(what2pl)
  if what2pl(pl) == 0 % ff, fr
    nexttile
    plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
    xlabel('time');
    ylabel('net return');
    hold on
    plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
    legend('ff','fr');
    hold off
  elseif what2pl(pl) == 1 %f, ff, fr, fb
    nexttile
    plot(t(plw(1):plw(2)),f(plw(1):plw(2)),'k');
    xlabel('time');
    ylabel('net return');
    hold on
    plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
    plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
    plot(t(plw(1):plw(2)),fb(plw(1):plw(2)),'b');
    legend('f','ff','fr','fb');
    hold off
  elseif what2pl(pl) == 2 % fr-ff & f
    nexttile
    plot(t(plw(1):plw(2)),(fr(plw(1):plw(2))-ff(plw(1):plw(2))));
    %plot(t(plw(1):plw(2)),f(plw(1):plw(2)),'k');
    xlabel('time');
    ylabel('net return');

  elseif what2pl(pl) == 3 % f
    nexttile 
    plot(plw(1):plw(2),f(plw(1):plw(2)),'k');
    legend('f');
  elseif what2pl(pl) == 4 % c
    %land cost
    nexttile
    plot(t(plw(1):plw(2)),C);%.2f\n',cor(1,2));
    legend('land cost');
  elseif what2pl(pl) == 5 % ff
    nexttile
    plot(t(plw(1):plw(2)),ff(plw(1):plw(2)),'g');
    legend('ff');
  elseif what2pl(pl) == 6 % fr
    nexttile
    plot(t(plw(1):plw(2)),fr(plw(1):plw(2)),'r');
    legend('fr');
  elseif what2pl(pl) == 7 % fb
    nexttile
    plot(t(plw(1):plw(2)),fb(plw(1):plw(2)),'r');
    legend('fb');
  elseif what2pl(pl) == 8 % fund
    nexttile 
    plot(t(plw(1):plw(2)),fundt);
    legend('fund');
    xlabel('time');
  elseif what2pl(pl) == 9 % tj over time
    nexttile
    plot(t(plw(1):plw(2)),tjs(plw(1):plw(2)))
    legend('tj')
    xlabel('time')
  elseif what2pl(pl) == 10 % tj histogram
    nexttile
    histogram(tjs);
    xlabel('tj');
  elseif what2pl(pl) == 11 % C histogram
    nexttile
    histogram(C(C~=0)); % used when plotting only distribution with certain t_j
    Ctitle = '';
    Ctitle = strcat(Ctitle, sprintf('af=%.1f,ar=%.1f, tj%d fr',param(8),param(15),cind));
    title(Ctitle)
    xlim([150,380]);
    ylim([0 620]);
    %histogram(C);
    xlabel('C');
  elseif what2pl(pl) == 12  % buying points
    nexttile
    scatter(t(plw(1):plw(2)),buy);
    legend('buying points')
  elseif what2pl(pl) == 13 % B histogram
    nexttile
    histogram(ben);
    xlabel('ben');
  else % mean of fr's 10 steps ahead
  	nexttile
  	frm = zeros(1,(plw(2)-plw(1)+1));
  	for i = plw(1):plw(2)
  		frm(i) = mean(fr(i:(i+10)));
  	end
  	plot(t(plw(1):plw(2)),frm(plw(1):plw(2)));
  	legend('fr fw mean');
    ylim([22 32]);
  	xlabel('time')
  end
end

foo = 1;