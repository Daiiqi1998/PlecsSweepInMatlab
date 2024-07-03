function [] = myBodeFix(gca_mag,gca_pha,f_MIN_bode,f_MAX_bode,FontSize,Option)
% Option = 'abs' or  'dB'
if nargin < 6
	Option = 'dB';
	if nargin < 5
		FontSize = 8;
	end
end


set(gca_mag,'FontName','Times New Roman','FontSize',FontSize,'color','none')
set(gca_mag,'GridAlpha',1,'MinorGridAlpha',1,'GridColor',[0.8,0.8,0.8],'MinorGridColor',[0.8,0.8,0.8]);
grid(gca_mag,'on')
set(gca_mag,"XScale","log",'XLim',[f_MIN_bode,f_MAX_bode],"XTickLabel","")
if strcmp(Option,'dB')
	ylabel(gca_mag,'\fontname{宋体}幅值\fontname{Times new roman}(dB)','Interpreter','tex','FontSize',FontSize)
else
	ylabel(gca_mag,'\fontname{宋体}幅值\fontname{Times new roman}(abs)','Interpreter','tex','FontSize',FontSize)
end

set(gca_pha,'FontName','Times New Roman','FontSize',FontSize,'color','none')
set(gca_pha,'GridAlpha',1,'MinorGridAlpha',1,'GridColor',[0.8,0.8,0.8],'MinorGridColor',[0.8,0.8,0.8]);
grid(gca_pha,'on')
set(gca_pha,"XScale","log",'XLim',[f_MIN_bode,f_MAX_bode])
% set(gca_pha,"YTick",[-180,-90,0,90,180])
ylabel(gca_pha,'\fontname{宋体}相位\fontname{Times new roman}(°)','Interpreter','tex','FontSize',FontSize)
xlabel('\fontname{宋体}频率\fontname{Times new roman}(Hz)','Interpreter','tex','FontSize',FontSize)
end

