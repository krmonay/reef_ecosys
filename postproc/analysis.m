% === ver 2017/03/09   Copyright (c) 2017 Takashi NAKAMURA  =====

evn_file = '.././output/eco5-env_his.csv';
eco_file = '.././output/eco5-ecosys_his.csv';

xmin=4; ymin=5;
PFDmax =2000;

env = readtable(evn_file,'Delimiter',',', 'ReadVariableNames', true);
eco = readtable(eco_file,'Delimiter',',', 'ReadVariableNames', true);

t_exp = (2+45/60)/24;  % day
t_stp = 3/24;      % day

tot_Pn = eco.coral1_Pn + eco.sedeco_Pn;
tot_G  = eco.coral1_G  + eco.sedeco_G ;
for i=1:8
    t_min = 4 + t_stp*(i-1);
    t_max = 4 + t_stp*(i-1) +t_exp;
    Pn(i) = mean(tot_Pn(t_min<eco.time & eco.time<t_max));
    G(i)  = mean(tot_G (t_min<eco.time & eco.time<t_max));
end