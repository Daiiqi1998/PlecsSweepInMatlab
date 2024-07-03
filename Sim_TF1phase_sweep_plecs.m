%% 在仿真中扫频测量某个端口的等效响应。
% plecs模型需存在同一路径，用4.7打开挂在后台。可以用变步长仿真。
% 在plecs模型中放一个signal outport,输出第一路为i(in),第二路为o(out)。
% 模型需要预设simulator parameters initialization中参数runT,injA,injF,
% runT设为模型运行时间，injA和injF设为小信号扰动源的幅值和频率。
% 需要电压电流扰动时，请活用受控源串并联，只支持同时设置一个小信号源。


%   注意每组输入输出之间的参考方向设置。这会影响到最后对Pha_T的计算。
% 2.在进行仿真时，可以预先设定好稳态初始值；
%   也可以等待进入稳态后再开始扫频。

%% 通用常量
path = pwd;     % 在当前路径定位plecs模型
model_name = 'Test1p';    % 模型名
save_name = 'Swp_A1';    % 扫频结果保存名。若缺省，则按当前时间保存
% 如果还需要调用别的路径下的自定义函数，请在这里 addpath

%% 仿真设置
% 计算需要扫到的频点范围。为了得到在log空间平均排布的整倍数，每段填入
% round(logspace(log10(最小频率),log10(最大频率),取点数量)/频率分辨率)*频率分辨率
% 比如：
F_FREQS = unique([...
    round( logspace(log10(5),   log10(40), 10)/1  )*1, ...
    round( logspace(log10(40),   log10(100), 15)/5  )*5, ...
    round( logspace(log10(100), log10(5e3), 30)/25  )*25]);
% 其中，unique表示不能重复；分频段按不同密度、不同频率分辨率取扫频点。优先取50Hz的整数倍，或整数分之一。

figure(1); semilogx(F_FREQS,ones(1,length(F_FREQS)),'*')   % 画个图验证扫频点排布，可以注释掉
F_FUND = 50;    % 基波频率
injA = 1;   % 注入小信号幅值
T_STEADY = 0.10;   % 认为从此开始，仿真模型处于稳态

k_t = F_FREQS/F_FUND;   % 计算扫频点相对基波的倍数
N_FREQS = length(k_t);   % 去掉重复的频点后，更新个数


%% 仿真
% 建线程，开模型
proxy = jsonrpc('http://localhost:1080', 'Timeout', 10);
proxy.plecs.load([path '/' model_name '.plecs']);
% proxy.plecs.scope([model_name '/Scope'], 'ClearTraces');

% 先扫一次获取稳态工作点矢量，以在小信号结果中减去它。
% 包括稳态参考模型在内，每个频点占用一次仿真（做不了一次性多频率注入）。
Comp_i = zeros(1,N_FREQS); Comp_o = zeros(1,N_FREQS);

% 测量稳态参考
K_PER = 1/F_FUND;   % 频点次数的分辨率
T_SCALE = round(1/K_PER)/F_FUND;  % 根据取频点的分辨率计算需要取的时间段长度

% 构造调用plecs时的仿真设置,struct里的参数取值会覆盖plecs模型内参数设置
simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', 0));
result = proxy.plecs.simulate(model_name, simStruct);   % 运行仿真

% 不管仿真是定步长、变步长，对采样线性插值就行
t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);    % 取稳态时间段

x = interp1(result.Time,result.Values(1,:),t);  
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_I = y(round(1+k_t/K_PER));     % 这个是换算完的FFT结果

x = interp1(result.Time,result.Values(2,:),t);
y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
Comp_O = y(round(1+k_t/K_PER));

% 再遍历每个频点扫频，测量小信号响应
simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', injA));
for ix = 1:length(F_FREQS)
    K_PER = gcd(F_FREQS(ix), F_FUND)/F_FUND;
    T_SCALE = round(1/K_PER)/F_FUND;  % 根据频点相对基波的分辨率，计算需要取的时间段长度。

    simStruct.ModelVars.injF = F_FREQS(ix);     % 对单次扫频再设置参数
    simStruct.ModelVars.runT = T_STEADY+T_SCALE;

    result = proxy.plecs.simulate(model_name, simStruct);   % 运行

    t = T_STEADY:0.1/F_FREQS(end):result.Time(end);     N_t = length(t);
    
    x = interp1(result.Time,result.Values(1,:),t);   
    y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
    Comp_i(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_I(ix);

    x = interp1(result.Time,result.Values(2,:),t);
    y=fft(x)*2/N_t; y=y(1:floor(N_t/2)+1);
    Comp_o(ix) = y(round(1+k_t(ix)/K_PER)) - Comp_O(ix);
end

Comp_T = Comp_o./Comp_i;
mag_T = abs(Comp_T);
pha_T = mod(phase(Comp_T)*180/pi + 180, 360)-180;  % 保证相位结果在(-180,+180)之内


%% 取名保存扫频结果
if isempty(save_name)
    tempt = datetime;
    save_name=string("SweepResult_"+...
	    yyyymmdd(tempt)+'_'+...
	    hour(tempt)+'_'+minute(tempt)+'_'+...
	    round(second(tempt)));
end
sweep_result=[F_FREQS;mag_T;pha_T];     % 注意幅值是以绝对值的形式给的，不是dB
eval("save('"+save_name+".mat', 'sweep_result');")

%% 画图
% 做阻抗的时候，我用绝对值画图。如果你希望用dB来画，那么请把这句改成true。
Plot_in_dB = false;

%{
% 如果已有扫频结果，就执行这一段来读取，不用重复扫了
F_FREQS = sweep_result(1,:);
mag_T = sweep_result(2,:);
pha_T = mod(sweep_result(3,:)+180,360)-180;
%}

f_MIN_bode = F_FREQS(1);
f_MAX_bode = F_FREQS(end);
%{
% 如果你有传函，可以把它定义在这里，或者事先读取它。比如：
s=tf('s'); 
sys_s = 1/(s*s+50*s+250);
[f_to_bode, mag_to_bode, pha_to_bode] = myBodeCal(f_MIN_bode, f_MAX_bode, 100, sys_s);
%}
TF_plot = false;



FontSize = 8;
gcf_bode=figure(2); clf; set(gcf,'Units','centimeters','Position',[2,8,8,6])

gca_mag=subplot(2,1,1); hold on
gca_pha=subplot(2,1,2); hold on


if Plot_in_dB
    plot_mag_1 = semilogx(gca_mag,F_FREQS,20*log10(mag_T),'k*');
    plot_pha_1 = semilogx(gca_pha,F_FREQS,pha_T,'k*');
    
    if TF_plot
        plot_mag_2 = semilogx(gca_mag,f_to_bode,mag_to_bode,'-','color',[0.8,0,0]);
        plot_pha_2 = semilogx(gca_pha,f_to_bode,pha_to_bode,'-','color',[0.8,0,0]);
        legend([plot_pha_1,plot_pha_2],'swp','model','FontSize',FontSize);	% 只选定需要标注图例的图像
    else
        legend([plot_pha_1],'swp','FontSize',FontSize);	% 只选定需要标注图例的图像
    end
    
    myBodeFix(gca_mag,gca_pha,f_MIN_bode,f_MAX_bode,FontSize,'dB')
    set(gca_mag,'YLim',[min(mag_to_bode),max(mag_to_bode)])
    set(gca_pha,'YLim',[min(pha_to_bode),max(pha_to_bode)])
else
    plot_mag_1 = semilogx(gca_mag,F_FREQS,mag_T,'k*');
    plot_pha_1 = semilogx(gca_pha,F_FREQS,pha_T,'k*');
    
    if TF_plot
        plot_mag_2 = semilogx(gca_mag,f_to_bode,10.^(mag_to_bode/20),'-','color',[0.8,0,0]);
        plot_pha_2 = semilogx(gca_pha,f_to_bode,pha_to_bode,'-','color',[0.8,0,0]);
        legend([plot_pha_1,plot_pha_2],'swp','model','FontSize',FontSize);	% 只选定需要标注图例的图像
    else
        legend([plot_pha_1],'swp','FontSize',FontSize);	% 只选定需要标注图例的图像
    end
    
    myBodeFix(gca_mag,gca_pha,f_MIN_bode,f_MAX_bode,FontSize,'abs')
    set(gca_mag,'YLim',[min(mag_T),max(mag_T)])
    set(gca_pha,'YLim',[min(pha_to_bode),max(pha_to_bode)])
end


%% 结束
%{
% 备份：matlab调用plecs仿真的例程里，是把所有仿真的数据都存下来再做处理的
    simStructs = cell(size(F_FREQS));
    % Initialize simStruct as cell array with all values for L1
    for ix = 1:length(F_FREQS)
       simStructs{ix}.ModelVars.injF = F_FREQS(ix);
    end
    results = proxy.plecs.simulate(model_name, simStructs);
    t = results(5).Time;
    x = results(5).Values;
%}

%{
% 备注：对FFT功能的说明
% 取的时间段需要包括首尾，即 N_t=length(t) 需要为奇数；
% t表达的时间长度取倒数即为频谱的分辨率（若t=0.02，则计算得y的相邻数据间跨度为50Hz）
% t的分辨率*2的倒数为频谱的频率上限。
    t=0:0.0001:1;   N_t = length(t);
    x = 2*sin(2*pi*50*t) + 1*cos(2*pi*150*t);
    y=fft(x)*2/N_t;     % 计算频谱
% 到这步，y对应的频点分别是[dc,1,2,...,50,...,5000,-5000,...,-50,...,1]
    y=y(1:round(N_t/2)+1);     % 截取单边频谱
    plot(0:length(y)-1,abs(y))      % 我恨matlab的索引为什么不是从0开始
%}

