% === ver 2017/03/07   Copyright (c) 2017 Takashi NAKAMURA  =====

box_file = '.././output/eco5-box_his.csv';
ch1_file = '.././output/eco5-crl1_his.csv';
ch2_file = '.././output/eco5-crl2_his.csv';
ca1_file = '.././output/eco5-crl1_ave.csv';
ca2_file = '.././output/eco5-crl2_ave.csv';

xmin=0; ymin=5;
PFDmax =2000;

d   = readtable(box_file,'Delimiter',',', 'ReadVariableNames', true);
ch1 = readtable(ch1_file,'Delimiter',',', 'ReadVariableNames', true);
ch2 = readtable(ch2_file,'Delimiter',',', 'ReadVariableNames', true);
ca1 = readtable(ca1_file,'Delimiter',',', 'ReadVariableNames', true);
ca2 = readtable(ca2_file,'Delimiter',',', 'ReadVariableNames', true);


figure;
subplot(4,2,1); 
plot(ch1.time, ch1.Pn,'b');
hold on
plot(ch1.time, ch1.G, 'r');
axis([xmin ymin  -0.3 0.3])
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,2);
plot(ch1.time, ch1.pHcal,'r');
axis([xmin ymin  7.5 9.5])
hold on
plot(ch1.time, ch1.pHcoe, 'g');
plot(ch1.time, ch1.pHamb, 'b');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,3);
plot(ch1.time, ch1.R,'b');
axis([xmin ymin  0 0.6])
hold on
plot(ch1.time, ch1.Pg, 'r');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,4);
plot(ch1.time, ch1.DOcoe, 'g');
axis([xmin ymin  0 600])
hold on
plot(ch1.time, ch1.DOamb, 'b');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,5);
plot(ch1.time, ch1.TAcal,'r');
axis([xmin ymin  1600 3200])
hold on
plot(ch1.time, ch1.TAcoe, 'g');
plot(ch1.time, ch1.TAamb, 'b');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,6);
plot(ch1.time, ch1.Wacal,'r');
axis([xmin ymin  0 20])
hold on
plot(ch1.time, ch1.Waamb, 'b');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,7);
plot(ch1.time, ch1.DICcal,'r');
axis([xmin ymin  1200 2200])
hold on
plot(ch1.time, ch1.DICcoe, 'g');
plot(ch1.time, ch1.DICamb, 'b');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,8);
plot(ch1.time, ch1.QC,'g');
axis([xmin ymin  0 25])
hold on
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])
