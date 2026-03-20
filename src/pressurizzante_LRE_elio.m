clear
close all
clc

%% ossidante e combustibile
M_tot = (1596.645)/2;
r = 4.2;

M_ox = (r/(r+1))*M_tot;
rho_ox = 1520;
M_f = (1/(r+1))*M_tot;
rho_f = 800;

V_ox = M_ox/rho_ox;
V_f = M_f/rho_f;

%check 
M_ox_check = 644.347345612134;
M_f_check = 153.416034669556;
err_perc_Mox = abs(M_ox_check-M_ox)/M_ox*100
err_perc_Mf = abs(M_f_check-M_f)/M_f*100

%% risoluzione

V_vect = [];
T_vect = [360:1:600];

for i=1:length(T_vect)

%% dati 
T_i = T_vect(i); %parametro, valore_minimo = T_f
p_f = 4.137e6;
T_f = 70+273.15; % T_tank valore intermedio 
R_he = 2078.6;
gamma = 1.66;

%% p iniziale
p_i = p_f*((T_i/T_f)^((gamma)/(gamma-1)));

%% pressurizzante
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

V_vect = [V_vect,V_he];
end

%% scelgo T
diff = V_vect(1)-V_vect(end); 
diff = 0.98*diff;          %variazione del 98% del volume di ingombro
ingombro = V_vect(1)-diff; %ingombro scelto
for i=1:length(V_vect)
    if (V_vect(i)-ingombro) < 0.001
        volume = V_vect(i)
        temperatura = T_vect(i)
        j = i;
        break
    end
end

%% dimensionamento del cilindro
H = 1; %impongo h=1m
A = V_vect(j)/H;
diametro = sqrt(4*A/pi)

%% grafico ingombro
figure
plot(T_vect,V_vect,'LineWidth',2)
grid on;
xlabel('temperatura')
ylabel('volume')
hold on
plot(T_vect(j),V_vect(j),'rx-','LineWidth',4)

