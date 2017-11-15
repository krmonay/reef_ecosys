% === Copyright (c) 2017 Takashi NAKAMURA  =====

zh1_file = '.././output/eco5-zoo1_his.csv';
% zh1_file = '.././output/eco5-zoo2_his.csv';

xmin=0; xmax=5;
PFDmax =2000;

zh1 = readtable(zh1_file,'Delimiter',',', 'ReadVariableNames', true);

%% 

figure('PaperSize',[20 30],...
    'OuterPosition',[0 0 1000 1050]);
%     'GraphicsSmoothing','off',...
%     'Color',[1 1 1],...

subplot(4,2,1); 
plot(zh1.time, zh1.Pg,'b');
hold on
plot(zh1.time, zh1.R, 'r');
axis([xmin xmax  -0.5e-3 1e-3])
ylabel('Pg, R (pmol cell^-^1 s^-^1)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('Pg','R','E', 'Location','southoutside','Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,2);
plot(zh1.time, zh1.QC,'b');
axis([xmin xmax  0 20])
hold on
ylabel('QC (pmol cell^-^1)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('QC','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,3);
plot(zh1.time, zh1.dens,'b');
axis([xmin xmax  0 1.5e6])
hold on
ylabel('Zoox. density (cell cm^-^2)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('Zoox. dens','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,4);
plot(zh1.time, zh1.DOcoe, 'g');
axis([xmin xmax  0 800])
hold on
plot(zh1.time, zh1.DOamb, 'b');
ylabel('DO (umol kg^-^1)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('DOcoe','DOamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,5);
plot(zh1.time, zh1.TAcal,'r');
axis([xmin xmax  1000 5000])
hold on
plot(zh1.time, zh1.TAcoe, 'g');
plot(zh1.time, zh1.TAamb, 'b');
ylabel('TA (umol kg^-^1)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('TAcal','TAcoe','TAamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,6);
plot(zh1.time, zh1.Wacal,'r');
axis([xmin xmax  0 20])
hold on
plot(zh1.time, zh1.Waamb, 'b');
ylabel('Omega arg')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('Wcal','Wcoe','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,7);
plot(zh1.time, zh1.DICcal,'r');
axis([xmin xmax  500 4000])
hold on
plot(zh1.time, zh1.DICcoe, 'g');
plot(zh1.time, zh1.DICamb, 'b');
ylabel('DIC (umol kg^-^1)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('DICcal','DICcoe','DICamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,8);
plot(zh1.time, zh1.QC,'g');
axis([xmin xmax  0 25])
hold on
ylabel('CH2O (umol cm^-^2)')
yyaxis right
plot(zh1.time, zh1.PFD, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin xmax  0 PFDmax])
legend('CH2O','E', 'Location','southoutside','Orientation','horizontal')
