%%
% Nume si prenume: TUDOR ADRIAN-RAUL-STEFANEL
%
clearvars
clc
%% Magic numbers (replace with received numbers)
m = 3;
n = 9;
%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300);
a2 = (1000+n*300);
b0 = (2.2+m+n)/5.5;
rng(m+10*n)
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)];
%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1; % simulation duration
gain = 15;
umin = -gain; umax = gain; % input saturation
ymin = -b0*gain/1.8; ymax = b0*gain/1.8; % output saturation
whtn_pow_in = 1e-9*5*(((m-1)*5+n)/5)/2; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(9); % input quantizer (DAC)
whtn_pow_out = 1e-8*5*(((m-1)*8+n)/5)/2; % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)
u_op_region = -(m+n/5)/2; % operating point
%% Input setup (can be changed/replaced/deleted)
u0 = 0;     % fixed
ust = 3;  % must be modified (saturation)
t1 = 12/a1; % recommended
%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2022b.slx");
t = out.tout;
u = out.u;
y = out.y;
plot(t,u,t,y)
shg


%% System identification
i1=5071;
i2=5955;
i3=10481;
i4=12010;

figure; plot(t,u, t,y); grid on; hold on; legend('u','y')
plot(t(i1), y(i1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t(i2), y(i2), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t(i3), y(i3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t(i4), y(i4), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);


 
u0=mean(u(i1:i2));
ust=mean(u(i3:i4));
y0=mean(y(i1:i2));
yst=mean(y(i3:i4))
 
K=(yst-y0)/(ust-u0)

%%
%%partea reala a polilor
i5 = 5681;
i6 = 10566;

figure; plot(t,u, t,y); grid on; hold on; legend('u','y')
plot(t(i5), y(i5), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t(i6), y(i6), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);


t_aux = t(i5:i6);
y_aux = abs(y(i5:i6) - yst);

figure
plot(t_aux,y_aux)

i7 = 336;
i8 = 1241;
i9 = 2080;

figure 
plot(t_aux,y_aux)
plot(t_aux(i7), y_aux(i7), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t_aux(i8), y_aux(i8), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t_aux(i9), y_aux(i9), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);

t_reg=t_aux([i7,i8,i9]);
y_reg=log(y_aux([i7,i8,i9]));
figure
plot(t_reg,y_reg)
A_reg = [sum(t_reg.^2),sum(t_reg);
         sum(t_reg), length(t_reg)];
B_reg = [sum(t_reg.*y_reg);sum(y_reg)];
theta = inv(A_reg)*B_reg;
Re=theta(1)
%%
%%partea imaginara a polilor
i11 = 6000;
i12 = 6850;
figure; plot(t,u, t,y); grid on; hold on; legend('u','y')
plot(t(i11), y(i11), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);
plot(t(i12), y(i12), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 10);

%%
%7460 inmultire fara 2
Tosc=2*(t(i12)-t(i11))
Im=2.*pi/Tosc
%0.0017
%%zeta, wn
wn=sqrt(Re^2+Im^2)
zeta= -Re/wn
%%Validare
A=[0,1;-wn^2,-2*zeta*wn];
B=[0; K*wn^2];
C=[1,0];
D=[0];
sys=ss(A,B,C,D)
ysim2=lsim(sys,u,t,[y(1),10000]);
figure;
plot(t,u,t,y,t,ysim2)
J=1/sqrt(length(t)*norm(y-ysim2))
eMPN=norm(y-ysim2)/norm(y-mean(y))*100