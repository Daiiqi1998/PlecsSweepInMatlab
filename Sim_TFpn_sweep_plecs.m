%% 在plecs仿真模型中对2*2序阻抗做扫频测量。这里是可以修改的副本。
% 模型需要预设参数runT,injA,injF（允许为负）,
% 其他要求同单相扫频一致。

%% 通用常量
path = pwd;     % 在当前路径定位plecs模型
model_name = 'TestPN';    % 模型名
save_name = 'SwpPN_Z0';    % 扫频结果保存名。若缺省，则按当前时间保存
% 如果还需要调用别的路径下的自定义函数，请在这里 addpath

%% 仿真设置
% 计算需要扫到的频点范围。仅取正半，不可含0。
F_FREQS = unique([...
    round( logspace(log10(5),   log10(40), 5)/1  )*1, ...
    round( logspace(log10(40),   log10(100), 5)/5  )*5, ...
    100 + round( logspace(log10(5),   log10(40), 5)/1  )*1, ...     % 这半是为了负序也能覆盖到低频
    100 + round( logspace(log10(40),   log10(100), 5)/5  )*5, ...
    round( logspace(log10(100), log10(5e3), 10)/25  )*25]);
F_FREQS = union(F_FREQS,2*50-F_FREQS);  % 所有需测的正负序输入输出频率。注意次同步正序的耦合输出仍是正序。

% figure(1); semilogx(F_FREQS,ones(1,length(F_FREQS)),'*')   % 画个图验证扫频点排布，可以注释掉
F_FUND = 50;    % 基波频率
injA = 1;   % 注入小信号幅值
T_STEADY = 1.00;   % 认为从此开始，仿真模型处于稳态
% 小问题：由于求耦合响应时需要看频率不同的双方的初始相位，所以T_STEADY需要为频率最小分辨率，即1Hz(1.0s)的整数倍。
% 否则，相位需要按照360°的分数倍进行额外加减。我懒得写了，所以就多仿一会儿吧。

k_t = F_FREQS/F_FUND;   % 计算扫频点相对基波的倍数
N_FREQS = length(k_t);   % 去掉重复的频点后，更新个数

%% 仿真
% 建线程，开模型
proxy = jsonrpc('http://localhost:1080', 'Timeout', 10);
proxy.plecs.load([path '/' model_name '.plecs']);
% proxy.plecs.scope([model_name '/Scope'], 'ClearTraces');

%% 无注入运行1次，获取稳态工作点矢量，后续需要在小信号结果中减去它。
% 考虑到谐波耦合的可能性，包括稳态参考模型在内，每个频点占用一次仿真（做不了一次性多频率注入）。
Comp_in = zeros(1,N_FREQS);
Comp_out_self = zeros(1,N_FREQS);
Comp_out_coup = zeros(1,N_FREQS);

% 测量稳态参考
K_PER = 1/F_FUND;   % 所有待测频点次数的预计最小分辨率。如果需要测1Hz以下的东西，要改它。
T_SCALE = round(1/K_PER)/F_FUND;  % 根据取频点的分辨率计算需要取的时间段长度

% 构造调用plecs时的仿真设置,struct里的参数取值会覆盖plecs模型内参数设置
simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', 0));
result = proxy.plecs.simulate(model_name, simStruct);   % 运行仿真

% 不管仿真是定步长、变步长，对采样线性插值就行
t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);    % 取稳态时间段

x = interp1(result.Time,result.Values(1,:),t);  
y=fft(x)*2/N_t;
temp_index = k_t;    % 预生成fft的读取位置
for i=1:length(k_t)
    if k_t(i)>=0
        temp_index(i) = round(1 + k_t(i)/K_PER);
    else
        temp_index(i) = round(1 + length(y) + k_t(i)/K_PER);
    end
end
Comp_IN = y(temp_index);     % 这个是换算完的FFT结果

x = interp1(result.Time,result.Values(2,:),t);
y=fft(x)*2/N_t;
Comp_OUT_A = y(temp_index);

x = interp1(result.Time,result.Values(3,:),t);
y=fft(x)*2/N_t;
Comp_OUT_B = y(temp_index);

x = interp1(result.Time,result.Values(4,:),t);
y=fft(x)*2/N_t;
Comp_OUT_C = y(temp_index);

a = exp(1j*2*pi/3); % 相序分解需要的旋转常量


%% 2次扫频：p轴注入，测量p, n轴响应；n轴注入，测量p, n轴响应。
% 再遍历每个频点扫频，测量小信号响应
simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', injA));
for ix = 1:N_FREQS
    if F_FREQS(ix) == 0
        % 不能注入直流。
        Comp_in(ix) = nan;
        Comp_out_self(ix) = nan;
        Comp_out_coup(ix) = nan;
        continue
    end
    K_PER = gcd(F_FREQS(ix), F_FUND)/F_FUND;
    T_SCALE = round(1/K_PER)/F_FUND;  % 根据频点相对基波的分辨率，计算需要取的时间段长度。
    % injF 和 100-injF 对应的分辨率是一样的。

    simStruct.ModelVars.injF = F_FREQS(ix);     % 对单次扫频再设置参数
    simStruct.ModelVars.runT = T_STEADY+T_SCALE;

    result = proxy.plecs.simulate(model_name, simStruct);   % 运行

    t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);
    
    % 输入提取
    x = interp1(result.Time,result.Values(1,:),t);
    y=fft(x)*2/N_t;
    if k_t(ix)>=0
        temp_index = round(1 + k_t(ix)/K_PER);
    else
        temp_index = round(1 + length(y) + k_t(ix)/K_PER);
    end
    Comp_in(ix) = y(temp_index) - Comp_IN(ix);
    
    % 自响应输出提取
    x = interp1(result.Time,result.Values(2,:),t);
    y=fft(x)*2/N_t;
    tempA = y(temp_index) - Comp_OUT_A(ix);

    x = interp1(result.Time,result.Values(3,:),t);
    y=fft(x)*2/N_t;
    tempB = y(temp_index) - Comp_OUT_B(ix);

    x = interp1(result.Time,result.Values(4,:),t);
    y=fft(x)*2/N_t;
    tempC = y(temp_index) - Comp_OUT_C(ix);

    Comp_out_self(ix) = (tempA + a*tempB + a*a*tempC)/3;    % 由于已经考虑负频率，故不再区分正负序
    if k_t(ix)<0
        Comp_out_self(ix) = conj(Comp_out_self(ix));     % 但要对负序相位抵消一个相位取负
    end

    % 耦合响应输出提取
    ix1 = find(abs( k_t - (2-k_t(ix)) )<1e-4);     % ix是存数据的位置，ix1是从信号里找响应的位置
    if k_t(ix1) == 1 || k_t(ix1) == 2
        Comp_out_coup(ix) = nan;   % 缺省数
    else
        x = interp1(result.Time,result.Values(2,:),t);
        y=fft(x)*2/N_t;
        if k_t(ix1)>=0
            temp_index = round(1 + k_t(ix1)/K_PER);
        else
            temp_index = round(1 + length(y) + k_t(ix1)/K_PER);
        end
        tempA = y(temp_index) - Comp_OUT_A(ix);

        x = interp1(result.Time,result.Values(3,:),t);
        y=fft(x)*2/N_t;
        tempB = y(temp_index) - Comp_OUT_B(ix);

        x = interp1(result.Time,result.Values(4,:),t);
        y=fft(x)*2/N_t;
        tempC = y(temp_index) - Comp_OUT_C(ix);

        Comp_out_coup(ix) = (tempA + a*tempB + a*a*tempC)/3;
        if k_t(ix)<0
            Comp_out_self(ix) = conj(Comp_out_self(ix));
        end
    end

end

Comp_T_self = Comp_out_self(~isnan(Comp_out_self))./Comp_in(~isnan(Comp_out_self));
mag_T_self = abs(Comp_T_self);
pha_T_self = mod(phase(Comp_T_self)*180/pi+180, 360)-180;

Comp_T_coup = Comp_out_coup(~isnan(Comp_out_coup))./Comp_in(~isnan(Comp_out_coup));
mag_T_coup = abs(Comp_T_coup);
pha_T_coup = mod(phase(Comp_T_coup)*180/pi+180, 360)-180;

% 按照实际对应的频率正负，拆分到相序里
F_FREQS_SELF = F_FREQS(~isnan(Comp_out_self));
temp_index = 0;
for i =1:length(F_FREQS_SELF)-1
    if sign(F_FREQS_SELF(i)) ~= sign(F_FREQS_SELF(i+1))
        temp_index = i;
    end
end
F_FREQS_PP = F_FREQS_SELF(temp_index+1 : end);
mag_T_pp = mag_T_self(temp_index+1 : end);
pha_T_pp = pha_T_self(temp_index+1 : end);
F_FREQS_NN = flip(-F_FREQS_SELF(1 : temp_index));
mag_T_nn = flip(mag_T_self(1 : temp_index));
pha_T_nn = flip(pha_T_self(1 : temp_index));

F_FREQS_COUP = F_FREQS(~isnan(Comp_out_coup));
temp_index = 0;
for i =1:length(F_FREQS_COUP)-1
    if sign(F_FREQS_COUP(i)) ~= sign(F_FREQS_COUP(i+1))
        temp_index = i;
    end
end
F_FREQS_PN = F_FREQS_COUP(temp_index+1 : end);
mag_T_pn = mag_T_coup(temp_index+1 : end);
pha_T_pn = pha_T_coup(temp_index+1 : end);
F_FREQS_NP = flip(-F_FREQS_COUP(1 : temp_index));
mag_T_np = flip(mag_T_coup(1 : temp_index));
pha_T_np = flip(pha_T_coup(1 : temp_index));


%% 取名保存扫频结果
if isempty(save_name)
    tempt = datetime;
    save_name=string("SweepResult_"+...
        yyyymmdd(tempt)+'_'+...
        hour(tempt)+'_'+minute(tempt)+'_'+...
        round(second(tempt)));
end
sweep_result={[F_FREQS_PP;mag_T_pp;pha_T_pp], ...
              [F_FREQS_PN;mag_T_pn;pha_T_pn], ...
              [F_FREQS_NP;mag_T_np;pha_T_np], ...
              [F_FREQS_NN;mag_T_nn;pha_T_nn]};
% 注意幅值是以绝对值的形式给的，不是dB
eval("save('"+save_name+".mat', 'sweep_result');")

%% 画图
%{
% 如果已有扫频结果，就执行这一段来读取，不用重复扫了
F_FREQS_PP = sweep_result{1}(1,:); mag_T_pp = sweep_result{1}(2,:); pha_T_pp = sweep_result{1}(3,:);
F_FREQS_PN = sweep_result{2}(1,:); mag_T_pn = sweep_result{2}(2,:); pha_T_pn = sweep_result{2}(3,:);
F_FREQS_NP = sweep_result{3}(1,:); mag_T_np = sweep_result{3}(2,:); pha_T_np = sweep_result{3}(3,:);
F_FREQS_NN = sweep_result{4}(1,:); mag_T_nn = sweep_result{4}(2,:); pha_T_nn = sweep_result{4}(3,:);
%}

f_MIN_bode = F_FREQS_PP(1); f_MAX_bode = F_FREQS_PP(end);

FontSize = 8;
gcf_bode=figure(2); clf; set(gcf,'Units','centimeters','Position',[2,6,16,12])

gca_mag_pp = subplot('Position', [0.06, 0.75, 0.42, 0.22]); hold on
gca_mag_pn = subplot('Position', [0.56, 0.75, 0.42, 0.22]); hold on
gca_mag_np = subplot('Position', [0.06, 0.27, 0.42, 0.22]); hold on
gca_mag_nn = subplot('Position', [0.56, 0.27, 0.42, 0.22]); hold on
gca_pha_pp = subplot('Position', [0.06, 0.53, 0.42, 0.22]); hold on
gca_pha_pn = subplot('Position', [0.56, 0.53, 0.42, 0.22]); hold on
gca_pha_np = subplot('Position', [0.06, 0.05, 0.42, 0.22]); hold on
gca_pha_nn = subplot('Position', [0.56, 0.05, 0.42, 0.22]); hold on

plot_mag_pp_1 = semilogx(gca_mag_pp,F_FREQS_PP,mag_T_pp,'-k*');
plot_pha_pp_1 = semilogx(gca_pha_pp,F_FREQS_PP,pha_T_pp,'-k*');
plot_mag_pn_1 = semilogx(gca_mag_pn,F_FREQS_PN,mag_T_pn,'-k*');
plot_pha_pn_1 = semilogx(gca_pha_pn,F_FREQS_PN,pha_T_pn,'-k*');
plot_mag_np_1 = semilogx(gca_mag_np,F_FREQS_NP,mag_T_np,'-k*');
plot_pha_np_1 = semilogx(gca_pha_np,F_FREQS_NP,pha_T_np,'-k*');
plot_mag_nn_1 = semilogx(gca_mag_nn,F_FREQS_NN,mag_T_nn,'-k*');
plot_pha_nn_1 = semilogx(gca_pha_nn,F_FREQS_NN,pha_T_nn,'-k*');

% legend([plot_pha_1],'swp','FontSize',FontSize);   % 只选定需要标注图例的图像

myBodeFix(gca_mag_pp,gca_pha_pp,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_pn,gca_pha_pn,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_np,gca_pha_np,f_MIN_bode,f_MAX_bode,FontSize,'abs')
myBodeFix(gca_mag_nn,gca_pha_nn,f_MIN_bode,f_MAX_bode,FontSize,'abs')
set(gca_mag_pp,'YLim',[min(mag_T_pp),max(mag_T_pp)])
set(gca_pha_pp,'YLim',[min(pha_T_pp),max(pha_T_pp)])
set(gca_mag_pn,'YLim',[min(mag_T_pn),max(mag_T_pn)])
set(gca_pha_pn,'YLim',[min(pha_T_pn),max(pha_T_pn)])
set(gca_mag_np,'YLim',[min(mag_T_np),max(mag_T_np)])
set(gca_pha_np,'YLim',[min(pha_T_np),max(pha_T_np)])
set(gca_mag_nn,'YLim',[min(mag_T_nn),max(mag_T_nn)])
set(gca_pha_nn,'YLim',[min(pha_T_nn),max(pha_T_nn)])

