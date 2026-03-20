clear 
close all
clc

%% dati
L = 2; %lunghezza caretteristica assunta con considerazioni dalle tabulate
d_c = 0.0883;
A_c = pi*(d_c/2)^2;
d_g = 0.0516;
A_g = pi*(d_g/2)^2;
p_c = 6894.76*525;
rho_ox = 1520;
rho_f = 800;
m_ox = 5.5788;
m_f = 1.3283;

%% camera di combustione
V_c = L*A_g;
L_c = V_c/A_c  
A_c_lat = pi*d_c*L_c;

%% testata di iniezione
C_D = 0.7; 
d_h = 0.0005; %short-tube with rounded entrance
delta_p = 0.1*p_c;

%area iniettori
A_ox = m_ox/(C_D*sqrt(2*rho_ox*delta_p));
A_f = m_f/(C_D*sqrt(2*rho_f*delta_p));

%numero iniettori
A_h = pi*(d_h/2)^2;
N_ox = ceil(A_ox/A_h);
N_f = ceil(A_f/A_h); 
N_h = N_f %scegliamo doppiette quindi n minore

%modifica aree
A_ox_new = A_ox/N_h;
d_ox_new = sqrt(4*A_ox_new/pi)
A_f_new = A_f/N_h;
d_f_new = sqrt(4*A_f_new/pi)

%angolo di iniezione
v_ox = C_D*sqrt(2*delta_p/rho_ox);
v_f = C_D*sqrt(2*delta_p/rho_f);

alpha_ox_vect = [1:1:30]; 
alpha_f_vect = [];
for i=1:length(alpha_ox_vect)
    alpha_ox = alpha_ox_vect(i)/180*pi;
    if (m_ox*v_ox*sin(alpha_ox))/(m_f*v_f) < 1
        alpha_f = asin((m_ox*v_ox*sin(alpha_ox))/(m_f*v_f));
        alpha_f = alpha_f/pi*180;
        alpha_f_vect = [alpha_f_vect,alpha_f];
    end
end

alpha_ox_vect(15) %scelgo 15 gradi
alpha_f_vect(15)


