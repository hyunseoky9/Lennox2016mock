function [cumb,fund,buy,notenough,crinotmet] = buystrat(buy,code,cumb,fund,don,B,c,E,xf,xr,al,be,cvalth,mc,Lxf,Hxf,Lxr,Hxr,Lx,Hx,me)
%% buying strategy (bstrat):
%% input
%% code
%% 1 = get cval (B/c)^alpha*(donation)^beta and buy above some threshold of cval
%% 2 = buy when cost low
%% 3 = buy when cost high
%% 4 = buy when xf low
%% 5 = buy when xf high
%% 6 = buy when xr low
%% 7 = buy when xr high
%% 8 = buy when E low
%% 9 = buy when E high
%% don = this year's donation
%% output
%% cumb = cumulative benefit
%% fund = remaining fund
crinotmet = 0;
notenough = 0;
if c <= fund
	notenough = 1;
end
if code == 1
	cval = (B/c)^al*don^be; % conservation value

	%cvalm = cvalm + cval;
	%if cval >= mcval(2)
	%  mcval(2) = cval;
	%end
	%if cval <= mcval(1)
	%  mcval(1) = cval;
	%end
	if cval >= cvalth
	  if c <= fund % buy
		  cumb = cumb + B;
		  fund = fund - c;
		  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	  else
	  	buy = [buy 0];
	  end
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 2 % buy when c low 
    %fprintf("c=%.2f, Lc=%.2f, fund=%.2f\n",c,Lc,fund);
	if c <= mc 
	  if c <= fund % buy
		  cumb = cumb + B;
		  fund = fund - c;
		  buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('code=%d, bought, B=%.2f\n',code, B);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 3 % buy when cost high
	if c >= mc
	    if c <= fund % buy
		  cumb = cumb + B;
		  fund = fund - c;
		  buy = [buy 1];
		else
		  	buy = [buy 0];
		end
  	  %fprintf('code=%d, bought, B=%.2f\n',code, B);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 4 % buy when xf low
	if xf <= Lxf
	    if c <= fund % buy
			cumb = cumb + B;
			fund = fund - c;
			buy = [buy 1];
		else
			buy = [buy 0];
		end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
		crinotmet = 1;
		buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 5 % buy when xf high
	if xf >= Hxf 
	  if c <= fund % buy
		cumb = cumb + B;
		fund = fund - c;
		buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 6 % buy when xr low
	if xr <= Lxr
	  if c <= fund % buy
		  cumb = cumb + B;
		  fund = fund - c;
		  buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
		crinotmet = 1;
	    buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 7 % buy when xr high
	if xr >= Hxr 
	  if c <= fund % buy
	    cumb = cumb + B;
	    fund = fund - c;
	    buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 8 % buy when roi low
	if E <= me 
	  if c <= fund % buy
	    cumb = cumb + B;
	    fund = fund - c;
	    buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end	
else % buy when roi high
	if E >= me
	  if c <= fund % buy
	    cumb = cumb + B;
	    fund = fund - c;
	    buy = [buy 1];
	  else
	  	buy = [buy 0];
	  end
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  crinotmet = 1;
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
end