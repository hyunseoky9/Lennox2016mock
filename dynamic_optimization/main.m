parambundle = {timeline,a,x0r,xsig2,lf,af,xf0r,xfsig2,lr,ar,xr0r,xrsig2,...
lb,ab,xb0r,xbsig2,al,be,godsimnum,period,lag,A,burnin,cvalth,simtime,fund,...
cumb,bfn,rho,del,efmu,efsig2,ermu,ersig2,ch,b_def,t_jmethod,intrate,code};

stratstr = {'CVAL','Lc','Hc','Lxf','Hxf','Lxr','Hxr','LE','HE','HB','LB','Hx','Lx'};
pname = {'timeline','a','x0r','xsig2','lf','af','xf0r','xfsig2','lr','ar','xr0r','xrsig2',...
'lb','ab','xb0r','xbsig2','al','be','godsimnum','period','lag','A','burnin','cvalth','simtime','fund',...
'cumb','bfn','rho','del','efmu','efsig2','ermu','ersig2','ch','b_def','t_jmethod','intrate','code'};
paramset = paramsetmaker(parambundle,pname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting and simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iwantplot = 0;
holdonplot = 0; % when plotting things on a same plot
what2pl = [1,11]; % plot index
tjv = [1,3,5,7];

if iwantplot
  % mxn plot where m is number of things to plot and n is number of parameter set
  if any(what2pl == 12)
    m = length(what2pl) - 1 + length(tjv);
  else
    m = length(what2pl);
  end
  n = length(paramset);
  if holdonplot == 1
    hold on
  end
end

% initialize data files
interestp = interested_parameters(parambundle); %param interested in displaying 
filenum = wr(interestp,parambundle,pname,1,1,1,0,0,msg);

%% Simulation
for i = 1:length(paramset)
  % simulation
  param = paramset{i};
  receptacle = mainsim(param); % receptacle = {strcubm,t,x,xf,xr,xb,tjs,C,fundt,ben,buy}

  % simulation
  wr(interestp,parambundle,pname,receptacle,param,1,1,filenum); % write strategy result
  wr(interestp,parambundle,pname,receptacle,param,2,0,filenum); % write simulation data

  % print outs
  printout(pname,param,interestp,stratstr,receptacle);

  %% plotting 
  %sgtitle('param');
  if iwantplot
    if holdonplot == 1
      bleh = plotting(receptacle,param,11,2);
    else
      %pl(what2pl,tjv,length(param),receptacle,param,i);
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

end

if holdonplot== 1
  hold off
end
