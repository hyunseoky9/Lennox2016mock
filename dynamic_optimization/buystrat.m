function [cumb,fund,buy] = buystrat(buy,code,cumb,fund,B,c,cvalth,Lc,Hc,Lff,Hff,Lfr,Hfr,Lf,Hf)
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
%% output
%% cumb = cumulative benefit
%% fund = remaining fund
if code == 1
	cval = (B/c)^al*nfb(1)^be; % conservation value

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

else % buy when fr high
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
end