function[f_to_bode, mag_to_bode, pha_to_bode] = myBodeCal(f_MIN_bode, f_MAX_bode, N_bode, sys)

% 画图：不设置y轴的限制。需要的话自己补去。

if sign(f_MIN_bode)<0 && sign(f_MAX_bode)<0
	f_to_bode = -logspace(log10(-f_MAX_bode),log10(-f_MIN_bode),N_bode);  % 负频率轴的场合
elseif sign(f_MIN_bode)>0 && sign(f_MAX_bode)>0
	f_to_bode = logspace(log10(f_MIN_bode),log10(f_MAX_bode),N_bode);   % 正频率轴的场合
else
	fprintf("\n不行\n")
end

[mag_to_bode,pha_to_bode,wout] = bode(sys,2*pi*f_to_bode);	% 求解图线
mag_to_bode = mag_to_bode(:,:);
mag_to_bode = 20*log10(mag_to_bode);	% 可选：把abs化为dB
pha_to_bode = pha_to_bode(:,:);
for i=1:length(pha_to_bode)
	while pha_to_bode(i)>180
		pha_to_bode(i)=pha_to_bode(i)-360;
	end
	while pha_to_bode(i)<-180
		pha_to_bode(i)=pha_to_bode(i)+360;
	end
end
