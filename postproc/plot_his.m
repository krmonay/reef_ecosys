% === ver 2017/03/07   Copyright (c) 2017 Takashi NAKAMURA  =====

box_file  = 'eco5-box_his.csv';
ch1_file = 'eco5-crl1_his.csv';
ch2_file = 'eco5-crl2_his.csv';
ca1_file = 'eco5-crl1_ave.csv';
ca2_file = 'eco5-crl2_ave.csv';

xmin=0; ymin=5;
PFDmax =2000;

d   = readtable(box_file,'Delimiter',',', 'ReadVariableNames', true);
ch1 = readtable(ch1_file,'Delimiter',',', 'ReadVariableNames', true);
ch2 = readtable(ch2_file,'Delimiter',',', 'ReadVariableNames', true);
ca1 = readtable(ca1_file,'Delimiter',',', 'ReadVariableNames', true);
ca2 = readtable(ca2_file,'Delimiter',',', 'ReadVariableNames', true);


figure;
subplot(4,2,1); 
plot(d.time, d.coral1_Pn,'b');
hold on
plot(d.time, d.coral1_G, 'r');
axis([xmin ymin  -0.3 0.3])
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])

hold off
subplot(4,2,3);
plot(d.time, d.coral1_R,'b');
axis([xmin ymin  0 0.6])
hold on
plot(d.time, d.coral1_Pg, 'r');
yyaxis right
plot(d.time, d.PFDsurf, 'Color', [1 0.6 0]);
axis([xmin ymin  0 PFDmax])


