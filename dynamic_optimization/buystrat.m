function [cumb,fund,buy] = buystrat(buy,code,cumb,fund,don,B,c,E,ff,fr,al,be,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf,LE,HE)
%% buying strategy (bstrat):
%% input
%% code
%% 1 = get cval (B/c)^alpha*(donation)^beta and buy above some threshold of cval
%% 2 = buy when cost low
%% 3 = buy when cost high
%% 4 = buy when ff low
%% 5 = buy when ff high
%% 6 = buy when fr low
%% 7 = buy when fr high
%% 8 = buy when E low
%% 0 = buy when E high
%% don = this year's donation
%% output
%% cumb = cumulative benefit
%% fund = remaining fund
if code == 1
	cval = (B/c)^al*don^be; % conservation value

	%cvalm = cvalm + cval;
	%if cval >= mcval(2)
	%  mcval(2) = cval;
	%end
	%if cval <= mcval(1)
	%  mcval(1) = cval;
	%end
	if cval >= cvalth && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 2 % buy when c low 
    %fprintf("c=%.2f, Lc=%.2f, fund=%.2f\n",c,Lc,fund);
	if c <= Lc && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 3 % buy when cost high
	if c >= Hc && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 4 % buy when ff low
	if ff <= Lff && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 5 % buy when ff high
	if ff >= Hff && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 6 % buy when fr low
	if fr <= Lfr && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end

elseif code == 7 % buy when fr high
	if fr >= Hfr && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end
elseif code == 8 % buy when roi low
	if E <= LE && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end	
else % buy when roi high
	if E >= LE && c <= fund % buy
	  cumb = cumb + B;
	  fund = fund - c;
	  buy = [buy 1];
	  %fprintf('bought at t=%d, cval=%.2f\n',i,cval);
	  %fprintf('remaining fund=%.2f\n',fund);
	else
	  buy = [buy 0];
	  %fprintf('cost=%.2f',c);
	  %fprintf('remaining fund=%.2f\n',fund);
	end	
end