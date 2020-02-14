%clear all

timeline = 10000; %10000;
%general economy
a = 0.8;
x0r = 30;
%x0 = x0r*(1-a);
xsig2 = 10;
%sig = xsig2*(1-a^2);


%forestry
lf = 0.7; %lambda
af = 0;
xf0r = 25;
%xf0 = xf0r*(1-lf*af);
xfsig2 = 2;
%sigf = xfsig2*(1-lf^2*af^2);
%factorf = lf*(1-af)*x0r/(1-af*lf);

%housing
lr = 0.7;
ar = 1;
xr0r = 25;
%xr0 = xr0r*(1-lr*ar);
xrsig2 = 2;
%sigr = xrsig2*(1-lr^2*ar^2);

%donation
lb = 0.7;
ab = 0.9;
xb0r = 25;
%xb0 = xb0r*(1-lb*ab);
xbsig2 = 2;
%sigb = xbsig2*(1-lb^2*ab^2);

% buy strategy stuff
code = [2,3];
al = 1;
be = 1;

godsimnum = 1;

period = 10;
lag = 0;
A = 1; % amplitude

burnin = 201; % burn in first few values of net return as they have not converged yet

cvalth = 0; % buying threshold for buystrat code 1

simtime = 9500; % number of buying opportunities
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
pn = length(paramset); % number of parameter sets.
if any(what2pl == 12)
  m = length(what2pl) - 1 + length(tjv);
else
  m = length(what2pl);
end
n = length(paramset);


% strategy result data
interestp = [6,10,17]; %param interested in displaying 
filename = strcat('./data/st.csv');
fileID0 = fopen(filename,'w');
for i = 1:length(parambundle)
  l = sprintf('%.3f,',parambundle{i});
  fprintf(fileID0,'# %s: %s\n',pname{i},l(1:end-1));
end
header = sprintf('%s,',pname{interestp});
header = strcat(header,'cumb\n');
fprintf(fileID0,header);


for i = 1:length(paramset)
  % actual simulation
  param = paramset{i};
  receptacle = mainsim(param);
  strcumb = receptacle{1};

  %% write out data
  % strategy result
  l = sprintf('%.2f,',[param(interestp),receptacle{1}]);
  l = strcat(l(1:end-1),'\n');
  fprintf(fileID0,l);
  % simulation data
  fileparam = sprintf('%.2f,',param(interestp));
  filename = strcat('./data/',fileparam,'.csv');
  fileID = fopen(filename,'w');
  st = param(26); %simtime
  fprintf(fileID,'#');
  fprintf(fileID,'%.3f  ',param);
  fprintf(fileID,'\n');
  fprintf(fileID,'t,x,xf,xr,xb,tjs,C,fundt,ben,buy\n');
  fprintf(fileID,'%d,%.2f,%.2f,%.2f,%.2f,%d,%.2f,%.2f,%.2f,%d\n',...
    [receptacle{2}(1:st);receptacle{3}(1:st);receptacle{4}(1:st);receptacle{5}(1:st);receptacle{6}(1:st);receptacle{7};receptacle{8};receptacle{9};receptacle{10};receptacle{11}]); %t,x,xf,xr,xb,fundt,tjs,C,ben,buy
  fclose(fileID);

  % print outs
  fprintf('ar=%.2f, af=%.2f\n',param(10),param(6));
  for j = 1:length(strcumb)
    fprintf('str: %s \t mean cumb:%.2f\n',stratstr{param(17)},strcumb);
  end
  fprintf('\n\n');

  %% making data files

  %% plotting 
  %sgtitle('param');
  if holdonplot == 1
    bleh = plotting(receptacle,param,11,2);
  else
    for j = 1:length(what2pl)
      if what2pl(j) == 12 % if contingent c dist has to be plotted for multiple contingents
        for k = 1:length(tjv)
          subplot(m,n,(pn*k-(pn-1))+(i-1))
          bleh = plotting(receptacle,param,what2pl(j),tjv(k));
        end
      else
        subplot(m,n,(pn*j-(pn-1))+(i-1))
        bleh = plotting(receptacle,param,what2pl(j),0); %make it better
      end
    end
  end
end

fclose(fileID0);
if holdonplot== 1
  hold off
end


%param = [timeline,a,x0r,xsig2,lf,af,xf0r,xfsig2,lr,ar,xr0r,xrsig2,...
%lb,ab,xb0r,xbsig2,code(1),al(1),be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
%cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def];
%receptacle = mainsim(param);

% data from last parameter set.
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
