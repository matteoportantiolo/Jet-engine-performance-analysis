clear
close all
clc

%% dati di progetto
PORTATA_A=21.2;
THRUST=8830;
SFC=0.1428;
eta_b=0.95;
T4_tot=1058.15;
eta_t=0.795;
beta_c=3.5;
eta_c=0.78;
phi_in=0.0254*20;
M_in=0.771018355445866;%vin=870 kmh
gdr_c=1;
gdr_t=0.2;
%% dati
%camera
DELTAH=42e6;
T_vect=[];
for pi_b=0.90:0.01:0.98;
%turbina

for eta_tm=0.90:0.01:0.98;
%compressore

eta_cm=eta_tm;

%ugello
for eta_n=0.9:0.01:0.98;
%diffusore

for pi_d=0.90:0.01:0.98;

%aria
M_mol=28.95; %g/mol
gamma_a=1.4;
cp_a=1004;
R_a=cp_a-cp_a/gamma_a;
%gas combusti
gamma_gc=1.33;
cp_gc=1155;
R_gc=cp_gc-cp_gc/gamma_gc;
%ambiente
z=6000;
g=9.81;
Pamb0=101325;
Tamb0=288.15;
a=0.0065;
Tamb=Tamb0-a*z;
Pamb=Pamb0*(1-(a*z)/Tamb0)^(g/(R_a*a));
rho_amb=Pamb/(R_a*Tamb);
%% soluzione ciclo

a_in=sqrt(gamma_a*R_a*Tamb);
v_in=M_in*a_in;
portata_a=rho_amb*(pi*(phi_in^2)/4)*v_in;
T1_tot=Tamb*(1+((gamma_a-1)/2)*M_in^2);
P1_tot=Pamb*(1+((gamma_a-1)/2)*M_in^2)^(gamma_a/(gamma_a-1));
P2_tot=P1_tot*pi_d;
T2_tot=T1_tot;
P3_tot=beta_c*P2_tot;
T3_tot_id=T2_tot*(P2_tot/P3_tot)^((1-gamma_a)/gamma_a);
T3_tot=T2_tot+(T3_tot_id-T2_tot)/eta_c;
P4_tot=pi_b*P3_tot;
f=(cp_gc*T4_tot-cp_a*T3_tot)/(eta_b*DELTAH-cp_gc*T4_tot);
portata_f=f*portata_a;
T5_tot=T4_tot-(portata_a*cp_a*(T3_tot-T2_tot))/(eta_cm*eta_tm*(portata_a+portata_f)*cp_gc);
T5_tot_id=T4_tot+(T5_tot-T4_tot)/eta_t;
P5_tot=P4_tot*(T4_tot/T5_tot_id)^(gamma_gc/(1-gamma_gc));
T6_tot=T5_tot;
%ipotesi OE
P6=Pamb;
T6_id=T5_tot*(P5_tot/P6)^((1-gamma_gc)/gamma_gc);
T6=T5_tot-eta_n*(T5_tot-T6_id);
a_e=sqrt(gamma_gc*R_gc*T6);
v_e=sqrt(2*cp_gc*(T5_tot-T6));
M_e=v_e/a_e;
if M_e>1%ipotesi no OE
M_e=1;
T6=T5_tot*(2/(gamma_gc+1));
T6_id=T5_tot-(T5_tot-T6)/eta_n;
P6=P5_tot*(T5_tot/T6_id)^((1-gamma_gc)/gamma_gc);
a_e=sqrt(gamma_gc*R_gc*T6);
v_e=a_e;
end

%% prestazioni
rho_e=P6/(R_gc*T6);
A_e=(portata_a+portata_f)/(v_e*rho_e);
T=(portata_a+portata_f)*v_e-portata_a*v_in+(P6-Pamb)*A_e;
Isp_a=T/(g*portata_a);
TSFC=(portata_f/T)*3600;%kg/hN
P_av=portata_f*DELTAH;
P_p=T*v_in;
P_dis=0.5*(portata_a+portata_f)*(v_e-v_in)^2;
P_j=P_p+P_dis;
eta_th=P_j/P_av;
eta_p=P_p/P_j;
eta_g=eta_th*eta_p;
T_vect=[T_vect;T];
end
end
end
end