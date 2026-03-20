clear
close all
clc

T_i = 490;
p_f = 5.938e6;
T_f = 290;
R_he = 2078.6;
gamma = 1.66;

p_i = p_f*((T_i/T_f)^((gamma)/(gamma-1)))

V_ox = 5.24;
V_f = 4.02;

V_he = 0;
i = 0;
toll = 10;
M_vett = [];

while (toll > 0.02)
    i = i+1;
    M_he = (p_f/(R_he*T_f))*(V_ox+V_f+V_he);
    V_he = (M_he/p_i)*T_i*R_he;
    M_vett = [M_vett,M_he];
    if (i>1)
        toll = M_vett(i)-M_vett(i-1);
    end
end

M_he
V_he
rho = M_he/V_he %densità elio = 0.1784 (0 gradi,1 atm)
                %densità esercitazione = 13