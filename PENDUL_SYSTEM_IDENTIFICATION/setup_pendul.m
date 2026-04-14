%%
% Nume si prenume: TUDOR ADRIAN-RAUL-STEFANEL
%
clearvars
clc
%% Magic numbers (replace with received numbers)
m = 3;
n = 9;
%% Process data and experiment setup (fixed, do not modify)
Ts = 500e-6; % fundamental step size
u_star = 1.2+n*0.075;
delta = 0.125;
delta_spab = 0.075;
umin = -5; umax = 5; % input saturation
ymin = -100; ymax = 100; % output saturation
g = 9.81;
% pendulum parameters
M = 0.8-n/48;
l = 1.2-m/24;
b = 0.3+m/24;
% measurement
c1 = 180/pi;
c2 = 4+n/2;
% (theta0,omega0)
rng(m+10*n)
x0_slx = [(n+3)/50,(-1)^(n+1)*m/20];
% input white noise power and sampling time
whtn_pow_in = 1e-10*(Ts*1e4)/2;
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(13); % input quantizer (DAC)
% output white noise power and sampling time
whtn_pow_out = 1e-3*Ts;
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(13); % output quantizer (ADC)
meas_rep = round(7+n/2); % data acquisition hardware sampling limitation
%% Input setup (can be changed/replaced/deleted)
t1 = 0.3;
tr=0.5*2;%timpul de urcare aproximativ *2:3 datorita oscilatiilor
N=4;%la alegere
p= round(tr/N/Ts);
DeltaT = p*(2^N-1)*Ts*2;
[input_LUT_dSPACE,Tfin] = generate_input_signal(Ts,t1,DeltaT,N,p,u_star,delta,delta_spab);
%% Data acquisition (use t, u, y to perform system identification)
out = sim("pendul_R2022b.slx");
t = out.tout;
u = out.u;
y = out.y;
subplot(211)
plot(t,u)
subplot(212)
plot(t,y)
shg

%% System identification
i1=8463;
i2=25859;
i3=34008;
i4=51184;
N=11; %%dam zoom pe cel mai ingust esantion din graficul 1 si 3 si numaram cu Data Tips
%%date de identificare
t_id=t(i1:N:i2); %%Adaugarea N ul reprezinta decimarea
u_id=u(i1:N:i2);
u_id = u_id - mean(u_id); %%eliminare DC pt a corecta eroarea
y_id=y(i1:N:i2);
y_id = y_id - mean(y_id);
%date de validare
t_vd=t(i3:N:i4);
u_vd=u(i3:N:i4);
u_vd = u_vd - mean(u_vd);
y_vd=y(i3:N:i4);
y_vd = y_vd - mean(y_vd);
figure;
subplot(221)
plot(t_id, u_id)
subplot(223)
plot(t_id,y_id)
subplot(222)
plot(t_vd,u_vd)
subplot(224)
plot(t_vd,y_vd)
dat_id = iddata(y_id, u_id, t_id(2)-t_id(1));
dat_vd = iddata(y_vd, u_vd, t_vd(2)-t_vd(1));
%%
model_armax = armax(dat_id, [3,3,4,2])
figure, resid(model_armax, dat_vd)
figure, compare(model_armax, dat_vd)
armax = ss(model_armax)

armax_pem = pem(dat_vd, model_armax);
figure, resid(armax_pem, dat_vd);
figure, compare(armax_pem, dat_vd);
%%
model_oe = oe(dat_id, [3,3,1])
figure, resid(model_oe, dat_vd)
figure, compare(model_oe, dat_vd)
oe = ss(model_oe)
%%
model_ssest = ssest(dat_id, 2)
figure, resid(model_ssest, dat_vd)
figure, compare(model_ssest, dat_vd)


