clear
close all
clc

D1 = 17;
D1 = D1*2.54/100;
A1 = pi*D1^2/4;
%A0_vett = 0.09:0.01:0.2;
cp = 1004;
gam = 1.4;
R = 287;
T0 = 223.15;
P0 = 26500;
rho0 = P0/(R*T0);
v0_vett = 200:50:550;
v0_vett = v0_vett.*1.609./3.6;
m_a = [6.73015895073036, 6.93429902022430, 7.17823243325070, ...
    7.46792575019269, 7.81400069765478, 8.21734233040033, 8.68214908256387, 9.21311476878258];
A0_dati = m_a./(rho0.*v0_vett);
a0 = sqrt(gam*R*T0);
M0 = v0_vett./a0;
T0_t = T0.*(1+(gam-1)./2.*M0.^2);
P0_t = P0.*(1+(gam-1)./2.*M0.^2).^(gam/(gam-1));
% assumiamo flusso adiabatico isoentropico prima della presa
P1_t = P0_t;
T1_t = T0_t;

v1 = zeros(length(v0_vett),1);
T1 = zeros(length(v0_vett),1);
a1 = zeros(length(v0_vett),1);
M1 = zeros(length(v0_vett),1);
P1 = zeros(length(v0_vett),1);
rho1 = zeros(length(v0_vett),1);

for i = 1:length(v0_vett)
    v0 = v0_vett(i);
    A0 = A0_dati(i);
    fun = @(v1) log(A1/A0) + log(v1/v0) + ...
                cp/(gam*R)*log((gam*R*T0_t-gam*R*v1^2/(2*cp))/(gam*R*T0_t-gam*R*v0^2/(2*cp)));
    v_guess = 100;
    v1(i) = fzero(fun, v_guess);
    T1(i) = T0_t(i) - v1(i)^2/(2*cp);
    a1(i) = sqrt(gam*R*T1(i));
    M1(i) = v1(i)/a1(i);
    P1(i) = P0_t(i)/(1+(gam-1)/2*M1(i)^2)^(gam/(gam-1));
    rho1(i) = P1(i)/(R*T1(i));

end

%% presa d'aria
pi_d = 0.97;
P2_t = pi_d.*P1_t;
T2_t = T1_t;
M2 = 0.35;
P2 = P2_t./(1+(gam-1)/2.*M2.^2).^(gam/(gam-1));
T2 = T2_t./(1+(gam-1)/2.*M2.^2);
a2 = sqrt(gam*R.*T2);
v2 = M2.*a2;

