clear
close all
clc

%%  dati Motor 109-718
r = 4.2; 
T = 12.26e3;
t_b = 110; 
p_c = 6894.76*525;

rho_ox = 1520;
rho_f = 900;

gamma = 1.27; %assunta 
Gamma = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)));
Ru = 8314;

%aria
gamma_a = 1.4;
R_a = 287;
cp_a = 1004;
MM_aria = 28.96;
%quota
z = 9144;
p_amb0 = 101325;
T_amb0 = 288.15;
rho_amb0 = 1.225;
a=0.0065;
g = 9.81;
T_amb=T_amb0-a*z;
p_amb=p_amb0*(1-(a*z)/T_amb0)^(g/(R_a*a));
rho_amb=p_amb/(R_a*T_amb);

p_e = p_amb; %OE
g_0 = 9.81;

%% parametro radice di T_c/MM
v_e = 1775; %dato in ingresso invece che in uscita (intervallo 1737.36-1828.8 m/s)
          
par = v_e/sqrt((2*gamma)/(gamma-1)*Ru*(1-(p_e/p_c)^((gamma-1)/gamma)));

T_c_vett = [1800:100:2800];
MM_vett = [];
for i = 1:length(T_c_vett)
    MM_vett(i) = T_c_vett(i)/(par^2);
end

%plot
figure 
plot(T_c_vett,MM_vett,'LineWidth',2)
grid on
hold on 
plot(MM_aria*par^2,MM_aria,'rx-','LineWidth',4)
xlabel('Temperatura in camera di combustione [K]')
ylabel('Massa molecolare dei gas combusti [g/mol]')

for i = 1:length(T_c_vett)
    if MM_vett(i)<MM_aria %vincolo imposto su MM
        MM = MM_vett(i)   
        T_c = T_c_vett(i)
    end
end

%un solo valore rispetta il vincolo
R = Ru/MM;

%% prestazioni 
I_sp = v_e/g_0   
I_sp_x = 180;
errore_I_sp = abs(I_sp-I_sp_x)./I_sp_x*100

c_car = sqrt(R*T_c)/Gamma; %non conosco l'area di gola
c_T = (I_sp*g_0)/c_car;

%% portata 
m_p = T/(I_sp*g_0)  
m_p_x = 7;
errore_m_p = abs(m_p-m_p_x)./m_p_x*100

%% efflusso
T_e = T_c*((p_e/p_c)^((gamma-1)/gamma));
a_e = sqrt(gamma*R*T_e);
M_e = v_e/a_e

%% aree
E = (1/M_e)*(((2/(gamma+1))*(1+((gamma-1)/2)*M_e^2))^((gamma+1)/(2*(gamma-1))));
A_g = T/(p_c*c_T);
d_g = sqrt(4*A_g/pi)
A_e = E*A_g;
d_e = sqrt(4*A_e/pi)

%% camera combustione
M_c = 0.2; %imposto
a_c = sqrt(gamma*R*T_c);
v_c = M_c*a_c;
rho_c = p_c/(R*T_c);
A_c = m_p/(rho_c*v_c);
d_c = sqrt(4*A_c/pi)

%% ossidante e combustibile
m_ox = (r/(r+1))*m_p;
m_f = (1/(r+1))*m_p;

M_ox = m_ox*t_b*1.05; 
M_f = m_f*t_b*1.05;   
M_tot = M_f+M_ox

V_ox = M_ox/rho_ox;
V_f = M_f/rho_f;
V_tot = V_f+V_ox;

%% impulsi
rho_av = M_tot/V_tot;
I_v = I_sp*rho_av
I_tot = I_sp*M_tot*g_0

