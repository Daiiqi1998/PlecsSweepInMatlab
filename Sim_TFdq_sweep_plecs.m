%% 在plecs仿真模型中对同步系2*2响应做扫频测量。
% 模型需要预设参数runT,injA[2],injF,
% 其他要求同单相扫频一致。

%% 通用常量
path = pwd;     % 在当前路径定位plecs模型
model_name = 'TestDQ';    % 模型名
save_name = 'SwpDQ_Z0';    % 扫频结果保存名。若缺省，则按当前时间保存
% 如果还需要调用别的路径下的自定义函数，请在这里 addpath

%% 仿真设置
% 计算需要扫到的频点范围。为了得到在log空间平均排布的整倍数，每段填入
% round(logspace(log10(最小频率),log10(最大频率),取点数量)/频率分辨率)*频率分辨率
% 比如：
F_FREQS = unique([...
    round( logspace(log10(5),   log10(40), 5)/1  )*1, ...
    round( logspace(log10(40),   log10(100), 5)/5  )*5, ...
    round( logspace(log10(100), log10(5e3), 10)/25  )*25]);
% 其中，unique表示不能重复；分频段按不同密度、不同频率分辨率取扫频点。优先取50Hz的整数倍，或整数分之一。

% figure(1); semilogx(F_FREQS,ones(1,length(F_FREQS)),'*')   % 画个图验证扫频点排布，可以注释掉
F_FUND = 50;    % 基波频率
injOn = 1;   % 对开启注入的频道，注入小信号幅值
T_STEADY = 0.10;   % 认为从此开始，仿真模型处于稳态

k_t = F_FREQS/F_FUND;   % 计算扫频点相对基波的倍数
N_FREQS = length(k_t);   % 去掉重复的频点后，更新个数


%% 仿真
% 建线程，开模型
proxy = jsonrpc('http://localhost:1080', 'Timeout', 10);
proxy.plecs.load([path '/' model_name '.plecs']);
% proxy.plecs.scope([model_name '/Scope'], 'ClearTraces');

%% 无注入运行1次，获取稳态工作点矢量，后续需要在小信号结果中减去它。
% 考虑到谐波耦合的可能性，包括稳态参考模型在内，我让每个注入频点占用一次仿真。
Comp_in = zeros(1,N_FREQS);
Comp_out_dd = zeros(1,N_FREQS);
Comp_out_dq = zeros(1,N_FREQS);
Comp_out_qd = zeros(1,N_FREQS);
Comp_out_qq = zeros(1,N_FREQS);

% 测量稳态参考
K_PER = 1/F_FUND;   % 频点次数的分辨率
T_SCALE = round(1/K_PER)/F_FUND;  % 根据取频点的分辨率计算需要取的时间段长度

% 构造调用plecs时的仿真设置,struct里的参数取值会覆盖plecs模型内参数设置
simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', [0, 0]));
result = proxy.plecs.simulate(model_name, simStruct);   % 运行仿真

% 不管仿真是定步长、变步长，对采样线性插值就行
t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);    % 取稳态时间段

x = interp1(result.Time,result.Values(1,:),t);  
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_IN_D = y(round(1+k_t/K_PER));     % 这个是换算完的FFT结果

x = interp1(result.Time,result.Values(2,:),t);  
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_IN_Q = y(round(1+k_t/K_PER));     % 这个是换算完的FFT结果

x = interp1(result.Time,result.Values(3,:),t);
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_OUT_D = y(round(1+k_t/K_PER));

x = interp1(result.Time,result.Values(4,:),t);
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_OUT_Q = y(round(1+k_t/K_PER));

%% 2次扫频：d轴注入，测量d, q轴响应；q轴注入，测量d, q轴响应。
% 再遍历每个频点扫频，测量小信号响应
for batch = 1:2
    if batch == 1
        simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', [injOn, 0]));
    else
        simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', [0, injOn]));
    end
    for ix = 1:length(F_FREQS)
        K_PER = gcd(F_FREQS(ix), F_FUND)/F_FUND;
        T_SCALE = round(1/K_PER)/F_FUND;  % 根据频点相对基波的分辨率，计算需要取的时间段长度。
    
        simStruct.ModelVars.injF = F_FREQS(ix);     % 对单次扫频再设置参数
        simStruct.ModelVars.runT = T_STEADY+T_SCALE;
    
        result = proxy.plecs.simulate(model_name, simStruct);   % 运行
    
        t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);
        
        % 以下 result.Values(i,:) 分别对应 plecs 模型里标号为 i 的 out port
        if batch == 1
            % 从d轴注入谐波
            x = interp1(result.Time,result.Values(1,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_in(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_IN_D(ix);

            x = interp1(result.Time,result.Values(3,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_out_dd(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_OUT_D(ix);
            
            x = interp1(result.Time,result.Values(4,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_out_dq(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_OUT_Q(ix);
        else
            % 从q轴注入谐波
            x = interp1(result.Time,result.Values(2,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_in(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_IN_Q(ix);

            x = interp1(result.Time,result.Values(3,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_out_qd(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_OUT_D(ix);

            x = interp1(result.Time,result.Values(4,:),t);
            y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
            Comp_out_qq(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_OUT_Q(ix);
        end

    end
end
Comp_T_dd = Comp_out_dd./Comp_in;
mag_T_dd = abs(Comp_T_dd);
pha_T_dd = mod(phase(Comp_T_dd)*180/pi + 180, 360)-180;  % 保证相位结果在(-180,+180)之内
Comp_T_dq = Comp_out_dq./Comp_in;
mag_T_dq = abs(Comp_T_dq);
pha_T_dq = mod(phase(Comp_T_dq)*180/pi + 180, 360)-180;
Comp_T_qd = Comp_out_qd./Comp_in;
mag_T_qd = abs(Comp_T_qd);
pha_T_qd = mod(phase(Comp_T_qd)*180/pi + 180, 360)-180;
Comp_T_qq = Comp_out_qq./Comp_in;
mag_T_qq = abs(Comp_T_qq);
pha_T_qq = mod(phase(Comp_T_qq)*180/pi + 180, 360)-180;


%% 取名保存扫频结果
if isempty(save_name)
    tempt = datetime;
    save_name=string("SweepResult_"+...
	    yyyymmdd(tempt)+'_'+...
	    hour(tempt)+'_'+minute(tempt)+'_'+...
	    round(second(tempt)));
end
sweep_result=[F_FREQS; ...
    mag_T_dd;pha_T_dd; ...
    mag_T_dq;pha_T_dq; ...
    mag_T_qd;pha_T_qd; ...
    mag_T_qq;pha_T_qq];
% 注意幅值是以绝对值的形式给的，不是dB
eval("save('"+save_name+".mat', 'sweep_result');")

%% 画图

%{
% 如果已有扫频结果，就执行这一段来读取，不用重复扫了
F_FREQS = sweep_result(1,:);
mag_T_dd = sweep_result(2,:);
pha_T_dd = mod(sweep_result(3,:)+180,360)-180;
mag_T_dq = sweep_result(4,:);
pha_T_dq = mod(sweep_result(5,:)+180,360)-180;
mag_T_qd = sweep_result(6,:);
pha_T_qd = mod(sweep_result(7,:)+180,360)-180;
mag_T_qq = sweep_result(8,:);
pha_T_qq = mod(sweep_result(9,:)+180,360)-180;
%}

f_MIN_bode = F_FREQS(1); f_MAX_bode = F_FREQS(end);

FontSize = 8;
gcf_bode=figure(2); clf; set(gcf,'Units','centimeters','Position',[2,6,16,12])

% 这里用一下自定义画图指令
gca_mag_dd = subplot('Position', [0.06, 0.75, 0.42, 0.22]); hold on
gca_mag_dq = subplot('Position', [0.56, 0.75, 0.42, 0.22]); hold on
gca_mag_qd = subplot('Position', [0.06, 0.27, 0.42, 0.22]); hold on
gca_mag_qq = subplot('Position', [0.56, 0.27, 0.42, 0.22]); hold on
gca_pha_dd = subplot('Position', [0.06, 0.53, 0.42, 0.22]); hold on
gca_pha_dq = subplot('Position', [0.56, 0.53, 0.42, 0.22]); hold on
gca_pha_qd = subplot('Position', [0.06, 0.05, 0.42, 0.22]); hold on
gca_pha_qq = subplot('Position', [0.56, 0.05, 0.42, 0.22]); hold on

plot_mag_dd_1 = semilogx(gca_mag_dd,F_FREQS,mag_T_dd,'-k*');
plot_pha_dd_1 = semilogx(gca_pha_dd,F_FREQS,pha_T_dd,'-k*');
plot_mag_dq_1 = semilogx(gca_mag_dq,F_FREQS,mag_T_dq,'-k*');
plot_pha_dq_1 = semilogx(gca_pha_dq,F_FREQS,pha_T_dq,'-k*');
plot_mag_qd_1 = semilogx(gca_mag_qd,F_FREQS,mag_T_qd,'-k*');
plot_pha_qd_1 = semilogx(gca_pha_qd,F_FREQS,pha_T_qd,'-k*');
plot_mag_qq_1 = semilogx(gca_mag_qq,F_FREQS,mag_T_qq,'-k*');
plot_pha_qq_1 = semilogx(gca_pha_qq,F_FREQS,pha_T_qq,'-k*');

% legend([plot_pha_1],'swp','FontSize',FontSize);	% 只选定需要标注图例的图像

myBodeFix(gca_mag_dd,gca_pha_dd,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_dq,gca_pha_dq,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_qd,gca_pha_qd,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_qq,gca_pha_qq,f_MIN_bode,f_MAX_bode,FontSize,'abs')
set(gca_mag_dd,'YLim',[min(mag_T_dd),max(mag_T_dd)])
set(gca_pha_dd,'YLim',[min(pha_T_dd),max(pha_T_dd)])
set(gca_mag_dq,'YLim',[min(mag_T_dq),max(mag_T_dq)])
set(gca_pha_dq,'YLim',[min(pha_T_dq),max(pha_T_dq)])
set(gca_mag_qd,'YLim',[min(mag_T_qd),max(mag_T_qd)])
set(gca_pha_qd,'YLim',[min(pha_T_qd),max(pha_T_qd)])
set(gca_mag_qq,'YLim',[min(mag_T_qq),max(mag_T_qq)])
set(gca_pha_qq,'YLim',[min(pha_T_qq),max(pha_T_qq)])

