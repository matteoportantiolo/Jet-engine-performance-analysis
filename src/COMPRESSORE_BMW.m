clear
close all
clc
%% BMW0003
%I DUE CICLI FOR ANNIDATI SERVONO PER TROVARE L'ANGOLO ALPHA DI INGRESSO AL
%PRIMO STADIO DEL COMPRESSORE. 
%SONO FATTI SUL RAGGIO MEDIO SU CUI DEFINIAMO UN GRADO DI REAZIONE DI 0.5
%LE MATRICI SCORRONO GLI STADI SULLA COLONNA E TUTTI I POSSIBILI ALPHA SULLA RIGA

ALPHA=pi/8:pi/100000:pi/5;
%diametri e raggi
d_ext=21.65*0.0254;
d_int_in=d_ext*0.6;
hubtotip_out=0.725;
n_stadi=7;
d_c_int_out=d_ext*hubtotip_out;
m=((d_c_int_out-d_int_in)/2)/(n_stadi-1);
% N.B. 
%r_c vettore dei raggi che include tutti i raggi delle hub, tutti i raggi 
%medi e tutti i raggi delle tip dei vari stadi, utile per le tabelle in
%appendice
%per ottenere grafici più 'smooth' definire direttamente
%r_c=d_int_in/2:1/1000:d_ext/2;
r_c=d_int_in/2:m:d_c_int_out/2;
for k=1:n_stadi
    r_int=d_int_in/2+(k-1)*m;
    r_med=(r_int+(d_ext/2))/2;
    r_c=[r_c, r_med];
end
r_c=[r_c, d_ext/2];
r_m=(d_ext+d_int_in)/4;
hubtotip_in=d_int_in/d_ext;
%beta stadio decrescente
BETA_C=3.09;
BETA_stadio=BETA_C^(1/n_stadi);
BETA_st=zeros(n_stadi,1);
BETA_st(4,1)=BETA_stadio;
BETA_st(3,1)=BETA_stadio*1.01;
BETA_st(2,1)=BETA_stadio*1.02;
BETA_st(1,1)=BETA_stadio*1.03;
BETA_st(5,1)=BETA_stadio/1.01;
BETA_st(6,1)=BETA_stadio/1.02;
BETA_st(7,1)=BETA_stadio/1.03;
%inizializzazione matrici
alpha1=zeros(n_stadi,length(ALPHA));
alpha2=zeros(n_stadi,length(ALPHA));
U=zeros(n_stadi,1);
beta1=zeros(n_stadi,length(ALPHA));
beta2=zeros(n_stadi,length(ALPHA));
T1_tot=zeros(n_stadi,length(ALPHA));
P1_tot=zeros(n_stadi,length(ALPHA));
P2_tot=zeros(n_stadi,length(ALPHA));
P2_tot_reale=zeros(n_stadi,length(ALPHA));
P1_tot_reale=zeros(n_stadi,length(ALPHA));
T3_tot=zeros(n_stadi,length(ALPHA));
T3_tot_id=zeros(n_stadi,length(ALPHA));
C1_z=zeros(n_stadi,length(ALPHA));
C2_z=zeros(n_stadi,length(ALPHA));
W1_z=zeros(n_stadi,length(ALPHA));
W2_z=zeros(n_stadi,length(ALPHA));
C1_th=zeros(n_stadi,length(ALPHA));
C2_th=zeros(n_stadi,length(ALPHA));
W1_th=zeros(n_stadi,length(ALPHA));
W2_th=zeros(n_stadi,length(ALPHA));
rapporto_eulero=zeros(n_stadi,length(ALPHA));

R=0.5;
cp_a=1004;
gamma_a=1.4;
om=9500*2*pi/60;
perdita = 0.98;
eta_c=0.78;

for j=1:length(ALPHA)
    alpha1(1,j)=ALPHA(j);
    C1=155/cos(alpha1(1,j));
    %dati di temperature e pressioni ottenuti dal ciclo termodinamico
    %completo
    T1_tot(1,j)=253.572055801363;
    P1_tot(1,j)=41854.3456861369;
    P1_tot_reale(1,j)=41854.3456861369; 
    T3_tot_reale=377.219823622588;%uscita compressore
    % ricavare eta stadio 
    fun=@(x) eta_c-(BETA_C^((gamma_a-1)/gamma_a)-1)/(BETA_C^((gamma_a-1)/(gamma_a*x))-1);
    eta_st0=0.8;
    eta_st=fzero(fun,eta_st0);
    %altra funzione di controllo per eta stadio 
    f=@(y) T3_tot_reale/T1_tot(1,j)-BETA_C^((gamma_a-1)/(gamma_a*y));
    eta_check=fzero(f,eta_st0);
   
    for i=1:n_stadi      
        U(i)=r_m*om;
        C1_z(i,j)=C1*cos(alpha1(i,j));
        C1_th(i,j)=C1*sin(alpha1(i,j));
        C1_vect=[C1_z(i,j);C1_th(i,j)];
        W1_vect=C1_vect-[0;U(i)];
        W1_z(i,j)=W1_vect(1);
        W1_th(i,j)=W1_vect(2);
        beta1(i,j)=atan(W1_th(i,j)/W1_z(i,j));
        beta2(i,j)=atan(((-R+0.5)*2*norm([0;U(i)]))/C1_z(i,j)-tan(alpha1(i,j)));%uso grado di reazione
        C2_z(i,j)=C1_z(i,j);%ipotesi velocità assiale costante
        W2_z(i,j)=W1_z(i,j);
        W2_th(i,j)=W2_z(i,j)*tan(beta2(i,j));
        W2_vect=[W2_z(i,j);W2_th(i,j)];
        C2_vect=W2_vect+[0;U(i)];
        C2_th(i,j)=C2_vect(2);
        alpha2(i,j)=atan(C2_th(i,j)/C2_z(i,j));
        alpha1(i+1,j)=alpha1(i,j);%ipotesi stadio ripetuto
        rapporto_eulero(i,j)=(1+(om*r_m*(C2_z(i,j)*tan(alpha2(i,j))-C1_z(i,j)*tan(alpha1(i,j))))/(cp_a*T1_tot(i,j)))^(gamma_a/(gamma_a-1));
        P2_tot(i,j)=rapporto_eulero(i,j)*P1_tot(i,j);
        P2_tot_reale(i,j)=(BETA_st(i)/perdita)*P1_tot_reale(i,j);
        T3_tot_id(i,j)=T1_tot(i,j)*(BETA_st(i)^((gamma_a-1)/gamma_a));
        T3_tot(i,j)=T1_tot(i,j)+((T3_tot_id(i,j)-T1_tot(i,j))/eta_st);
        T1_tot(i+1,j)=T3_tot(i,j);
        P1_tot(i+1,j)=P2_tot(i,j);
        P1_tot_reale(i+1,j)=P2_tot_reale(i,j)*perdita; 
    end
end
P3_tot_reale=P1_tot_reale(2:end,:);
alpha3=alpha1(2:end,:);
%taglio l'ultima riga delle matrici '1'
P1_tot_reale=P1_tot_reale(1:end-1,:);
T1_tot=T1_tot(1:end-1,:);
alpha1=alpha1(1:end-1,:);
% IL CRITERIO PER SCEGLIERE ALPHA SELEZIONA L'ALPHA IN INGRESSO 
% ASSOCIATO AL RAPPORTO DI EULERO SUPERIORE MA PIU' VICINO
% AL BETA STADIO/PERDITA.
errore=rapporto_eulero-(BETA_st/perdita).*ones(1,length(ALPHA));
err=zeros(n_stadi,1);
index=zeros(n_stadi,1);
for i=1:n_stadi
    Err=errore(i,:);    
    [err(i),index(i)]=min(Err(Err>0));
end
[index,index_err]=min(index);
err=err(index_err);
errore_temperature=100*abs(T3_tot(end,index)-T3_tot_reale)/T3_tot_reale;
%% 
% QUI ENTRA IN GIOCO L'IPOTESI DI FREE VORTEX, OVVERO R*C_THETA COSTANTE
% LUNGO LA PALETTA
K1=C1_th(1,index)*r_m;
K2=C2_th(1,index)*r_m;
k=(K1+K2)/(2*om);%può essere usato per trovare i gradi di reazione(appendice)
%inizializzazione matrici
%'_r' INDICA CHE QUESTE SONO LE MATRICI CHE SCORRONO LE PALETTE
%IL NUMERO DI RIGHE RAPPRESENTA I RAGGI, MENTRE LE COLONNE GLI
%STADI
alpha1_r=zeros(length(r_c),n_stadi);
alpha2_r=zeros(length(r_c),n_stadi);
U_r=zeros(length(r_c),1);
beta1_r=zeros(length(r_c),n_stadi);
beta2_r=zeros(length(r_c),n_stadi);
c1z=155;
c2z=155;
c1_z_r=c1z*ones(length(r_c),n_stadi);
c2_z_r=c2z*ones(length(r_c),n_stadi);
w1_z_r=zeros(length(r_c),n_stadi);
w2_z_r=zeros(length(r_c),n_stadi);
c1_th_r=zeros(length(r_c),n_stadi);
c2_th_r=zeros(length(r_c),n_stadi);
w1_th_r=zeros(length(r_c),n_stadi);
w2_th_r=zeros(length(r_c),n_stadi);
grado_reazione_r=zeros(length(r_c),n_stadi); %MATRICE DEI GRADI DI REAZIONE
rapporto_eulero_r=zeros(length(r_c),n_stadi); %MATRICE DEI RAPPORTI ROTORE

for j=1:n_stadi
    for i=1:length(r_c)
        r=r_c(i);
        U_r(i)=om*r;
        c2_th_r(i,j)=K2/r;
        c2_vect=[c2z;c2_th_r(i,j)];
        w2_vect=c2_vect-[0;U_r(i)];
        w2_z_r(i,j)=w2_vect(1);
        w2_th_r(i,j)=w2_vect(2);
        beta2_r(i,j)=atan(w2_th_r(i,j)/w2_z_r(i,j));
        alpha2_r(i,j)=atan(c2_th_r(i,j)/c2z);
        c1_th_r(i,j)=K1/r;
        alpha1_r(i,j)=atan(c1_th_r(i,j)/c1z);
        c1_vect=[c1z;c1_th_r(i,j)];
        w1_vect=c1_vect-[0;U_r(i)];
        w1_z_r(i,j)=w1_vect(1);
        w1_th_r(i,j)=w1_vect(2);
        beta1_r(i,j)=atan(w1_th_r(i,j)/w1_z_r(i,j));
        grado_reazione_r(i,j)=0.5-c1z/U_r(i)*(tan(alpha1_r(i,j))+tan(beta2_r(i,j)))/2;
        rapporto_eulero_r(i,j)=(1+(om*r*(c2z*tan(alpha2_r(i,j))- ...
            c1z*tan(alpha1_r(i,j))))/(cp_a*T1_tot(j,1)))^(gamma_a/(gamma_a-1));
    end
end
%% 
%CONTROLLO MACH INGRESSO ROTORE
T1_tot_r=ones(length(r_c),n_stadi);
c1_r=zeros(length(r_c),n_stadi);
w1_r=zeros(length(r_c),n_stadi);
for i=1:n_stadi
    for j=1:length (r_c)
        T1_tot_r(j,i)=T1_tot(i,1);
    end
end
for i=1:n_stadi
    for j=1:length (r_c)
        c1_r(j,i)=sqrt((c1_th_r(j,i)).^2+(c1_z_r(j,i)).^2);
    end
end
for i=1:n_stadi
    for j=1:length (r_c)
        w1_r(j,i)=sqrt((w1_th_r(j,i)).^2+(w1_z_r(j,i)).^2);
    end
end
T1_s_r=T1_tot_r-c1_r.^2/(2*cp_a);%statica
a_r=sqrt(gamma_a*286.6.*T1_s_r);
Mach_r=c1_r./a_r;%assoluto
Mach_r_r=w1_r./a_r;%relativo

%% GRAFICI----utilizzare r_c con passo più fitto per aumentare la qualità
%tornare alla riga 24, eliminare da riga 24 a riga 30 e ridefinire r_c

%% PLOTTIAMO L'ANDAMENTO DELLE PRESSIONI TOTALI LUNGO I VARI STADI
tt = 0:0.5:n_stadi;
pressioni = [];
for i = 1:n_stadi
    pressioni = [pressioni, P1_tot(i,index)];
    pressioni = [pressioni, P2_tot(i,index)];
end
pressioni = [pressioni, P1_tot(end, index)];
pressioni_real = []; 
for i = 1:n_stadi
    pressioni_real = [pressioni_real, P1_tot_reale(i,index)];
    pressioni_real = [pressioni_real, P2_tot_reale(i,index)];
end
pressioni_real = [pressioni_real, P3_tot_reale(end, index)];
figure
plot(tt, pressioni,'linewidth',2)
hold on
plot(tt, pressioni_real,'linewidth',2)
xlabel('Stadio del compressore')
ylabel('Pressione totale [Pa]')
title('BMW 003')
legend('Eulero','Reale')

%% PLOTTIAMO LA VARIAZIONE DEGLI ANGOLI IN FUNZIONE DI r/rtip 
figure;
rateo=r_c./r_c(end);
plot(rateo,alpha1_r(:,1).*180/pi,'linewidth',2);
hold on
plot(rateo,beta1_r(:,1).*180/pi,'linewidth',2);
hold on 
plot(rateo,alpha2_r(:,1).*180/pi,'linewidth',2);
hold on
plot(rateo,beta2_r(:,1).*180/pi,'linewidth',2);
legend('alpha1','beta1','alpha2','beta2');
xlabel('r/r_t_i_p')
ylabel('angolo [°]')
title ('BMW 003')
%% PLOTTIAMO IL GRADO DI REAZIONE IN FUNZIONE DI r/rtip (AUMENTA DALL'HUB AL TIP)
figure;
plot(rateo,grado_reazione_r(:,1),'linewidth',2)
xlabel('r/r_t_i_p')
ylabel('R')
title('BMW 003')
%% PLOTTIAMO IL RAPPORTO PT2/PT1 LUNGO GLI STADI (DIMINUISCE)
figure;
X=1:0.01:n_stadi;
x=1:1:n_stadi;
Y=spline(x,rapporto_eulero_r(1,:),X);
plot(X,Y,'linewidth',2);
xlabel('Stadio del compressore')
ylabel('P_2_t_o_t/P_1_t_o_t')
title('BMW 003')
%% PLOTTIAMO IL RAPPORTO PT2/PT1 IN FUNZIONE DI ALPHA1
figure
plot(ALPHA,rapporto_eulero(1,:),'linewidth',2);
hold on
plot(alpha1(1,index),rapporto_eulero(1,index),'rx-','linewidth',4)
xlabel('alpha_1')
ylabel('P_2_t_o_t/P_1_t_o_t')
title('BMW 003 raggio medio primo stadio')
%% PLOTTIAMO IL MACH RELATIVO LUNGO I RAGGI
figure
plot(rateo, Mach_r_r(:,1),'linewidth',2);
hold on
plot(rateo, Mach_r_r(:,2),'linewidth',2);
plot(rateo, Mach_r_r(:,3),'linewidth',2);
plot(rateo, Mach_r_r(:,4),'linewidth',2);
plot(rateo, Mach_r_r(:,5),'linewidth',2);
plot(rateo, Mach_r_r(:,6),'linewidth',2);
plot(rateo, Mach_r_r(:,7),'linewidth',2);
legend('stadio 1','stadio 2','stadio 3','stadio 4','stadio 5','stadio 6','stadio 7')
xlabel('r/r_t_i_p')
ylabel('Ma_r_e_l')
title('BMW 003')