%clear all

timeline = 800; %10000;
%general economy
a = 0.8;
x0r = 30;
xsig2 = 10;

%forestry
lf = 0.7; %lambda
af = linspace(0,1,11);
xf0r = 25;
xfsig2 = 2;

%housing
lr = 0.7;
ar = linspace(0,1,11);
xr0r = 25;
xrsig2 = 2;

%donation
lb = 0.7;
ab = 0.9;
xb0r = 25;
xbsig2 = 2;

% buy strategy stuff
code = [1,2,3];
al = 1;
be = 1;

godsimnum = 1000;

period = 10;
lag = 0;
A = 1; % amplitude

burnin = 201; % burn in first few values of net return as they have not converged yet

cvalth = 0; % buying threshold for buystrat code 1

simtime = 300; % number of buying opportunities
fund = 0; % money saved
cumb = 0; % cummulative conservation value

bfn = 1; % benefit fn scheme. 1=constant, 2=normal var correlated to e_fj, 3=non-linear fn of 2

rho = 0.1; %0.1; % economic discount rate
del = 0.01; % ecological discount rate

efmu = 0; % mean of ef
efsig2 = 0; % var of ef
ermu = 0; % mean of er
ersig2 = 0; % var of er

ch = 0; % land change cost and option value

b_def = 10; % default b

parambundle = {timeline,a,x0r,xsig2,lf,af,xf0r,xfsig2,lr,ar,xr0r,xrsig2,...
lb,ab,xb0r,xbsig2,code,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def};
paramset = paramsetmaker(parambundle);

stratstr = {'CVAL','Lc','Hc','Lxf','Hxf','Lxr','Hxr','LE','HE'};
pname = {'timeline','a','x0r','xsig2','lf','af','xf0r','xfsig2','lr','ar','xr0r','xrsig2',...
'lb','ab','xb0r','xbsig2','code','al','be','godsimnum','period','lag','A','burnin','cvalth','simtime','fund',...
'cumb','bfn','rho','del','efmu','efsig2','ermu','ersig2','ch','b_def'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting and simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

holdonplot = 0; % when plotting things on a same plot
if holdonplot == 1
  hold on
end
what2pl = [12]; % plot index
tjv = [1,3,5,7];
if any(what2pl == 12)
  m = length(what2pl) - 1 + length(tjv);
else
  m = length(what2pl);
end
n = length(paramset);

% strategy result data
interestp = [6,10,17]; %param interested in displaying 
wr(interestp,parambundle,pname,1,1,1,0);

for i = 1:length(paramset)
  % simulation
  param = paramset{i};
  receptacle = mainsim(param); % receptacle = {strcubm,t,x,xf,xr,xb,tjs,C,fundt,ben,buy}

  % simulation
  wr(interestp,parambundle,pname,receptacle,param,1,1); % write strategy result
  wr(interestp,parambundle,pname,receptacle,param,2,0); % write simulation data

  % print outs
  fprintf('ar=%.2f, af=%.2f, str= %s\n',param(10),param(6),stratstr{param(17)});
  for j = 1:length(strcumb)
    fprintf('mean cumb:%.2f\n',receptacle{1});
  end
  fprintf('\n\n');

  %% plotting 
  %sgtitle('param');
  if holdonplot == 1
    bleh = plotting(receptacle,param,11,2);
  else
    pl(what2pl,tjv,length(param),receptacle,param,i);
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
  end

end

if holdonplot== 1
  hold off
end
