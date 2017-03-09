% === ver 2017/03/09   Copyright (c) 2017 Takashi NAKAMURA  =====

evn_file = '.././output/site04-env_his.csv';
ch1_file = '.././output/site04-crl1_his.csv';
% ch2_file = '.././output/site06-crl2_his.csv';
% ca1_file = '.././output/eco5-crl1_ave.csv';
% ca2_file = '.././output/eco5-crl2_ave.csv';

xmin=4; ymin=5;
PFDmax =2000;

env = readtable(evn_file,'Delimiter',',', 'ReadVariableNames', true);
ch1 = readtable(ch1_file,'Delimiter',',', 'ReadVariableNames', true);
% ch2 = readtable(ch2_file,'Delimiter',',', 'ReadVariableNames', true);
% ca1 = readtable(ca1_file,'Delimiter',',', 'ReadVariableNames', true);
% ca2 = readtable(ca2_file,'Delimiter',',', 'ReadVariableNames', true);


figure('PaperSize',[20 30],...
    'OuterPosition',[0 0 1000 1050]);
%     'GraphicsSmoothing','off',...
%     'Color',[1 1 1],...

subplot(4,2,1); 
plot(ch1.time, ch1.Pn,'b');
hold on
plot(ch1.time, ch1.G, 'r');
axis([xmin ymin  -0.3 0.3])
ylabel('G, Pn (nmol cm^-^2 s^-^1)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('Pn','G','E', 'Location','southoutside','Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,2);
plot(ch1.time, ch1.pHcal,'r');
axis([xmin ymin  7.5 9.5])
hold on
plot(ch1.time, ch1.pHcoe, 'g');
plot(ch1.time, ch1.pHamb, 'b');
ylabel('pH')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('pHcal','pHcoe','pHamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,3);
plot(ch1.time, ch1.R,'b');
axis([xmin ymin  0 0.6])
hold on
plot(ch1.time, ch1.Pg, 'r');
ylabel('Pg, R (nmol cm^-^2 s^-^1)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('R','Pg','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,4);
plot(ch1.time, ch1.DOcoe, 'g');
axis([xmin ymin  0 600])
hold on
plot(ch1.time, ch1.DOamb, 'b');
ylabel('DO (umol kg^-^1)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('DOcoe','DOamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,5);
plot(ch1.time, ch1.TAcal,'r');
axis([xmin ymin  1600 3200])
hold on
plot(ch1.time, ch1.TAcoe, 'g');
plot(ch1.time, ch1.TAamb, 'b');
ylabel('TA (umol kg^-^1)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('TAcal','TAcoe','TAamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,6);
plot(ch1.time, ch1.Wacal,'r');
axis([xmin ymin  0 20])
hold on
plot(ch1.time, ch1.Waamb, 'b');
ylabel('Omega arg')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('Wcal','Wcoe','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,7);
plot(ch1.time, ch1.DICcal,'r');
axis([xmin ymin  1200 2200])
hold on
plot(ch1.time, ch1.DICcoe, 'g');
plot(ch1.time, ch1.DICamb, 'b');
ylabel('DIC (umol kg^-^1)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('DICcal','DICcoe','DICamb','E', 'Location','southoutside','Orientation','horizontal')

hold off
subplot(4,2,8);
plot(ch1.time, ch1.QC,'g');
axis([xmin ymin  0 25])
hold on
ylabel('CH2O (umol cm^-^2)')
yyaxis right
plot(env.time, env.PFDsurf, 'Color', [1 0.6 0]);
ax = gca; ax.YColor = 'k';
ylabel('E (umol m^-^2 s^-^1)')
axis([xmin ymin  0 PFDmax])
legend('CH2O','E', 'Location','southoutside','Orientation','horizontal')
