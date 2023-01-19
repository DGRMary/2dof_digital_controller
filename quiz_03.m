clear all 
close all 
clc 

Ts = 10*10^-3;
s = tf('s');
G_cont = 30*(s+4)/(s*(s^2+7.5*s+15));
G = zpk(c2d(G_cont,Ts,'zoh'))
[z_G,p_G,k_G] = zpkdata(G,'v')
B = k_G*conv([1 -z_G(1)],[1 -z_G(2)]);
A = conv([1 -p_G(2)],[1 -p_G(3)]);

s_hat = 0.12;
ts_1 = 0.15;

zeta = 1.5*abs(log(s_hat))/sqrt(pi^2+(log(s_hat))^2)
wn = 4.6/(zeta*ts_1)

figure(1) 
zgrid(zeta,Ts*wn,'new')
zplane(zero(G),pole(G))

B_plus = [1 -z_G(2)];
B_minus = k_G*[1 -z_G(1)];
A_plus = A;
A_minus = 1;

p1c = -wn*zeta+j*wn*sqrt(1-zeta^2)
p2c = -wn*zeta-j*wn*sqrt(1-zeta^2)
p3c = -5*wn*zeta

p1 = exp(Ts*p1c)
p2 = exp(Ts*p2c)
p3 = exp(Ts*p3c)

Am = poly([p1 p2 p3])

A_dioph = conv([1 -1],conv([1 -1],A_minus));
B_dioph = B_minus;

[R1,S1,Am_check] = dioph_mtx(A_dioph,B_dioph,Am);
Am_check
S = conv(A_plus,S1);
R = conv([1 -1],conv(B_plus,R1));
C = zpk(tf(S,R,Ts))

G = zpk(G);
L = minreal(C*G);
Wprime = zpk(minreal(L/(1+L),10^-3));
kw = 1;
T_tilde =[1 -p3];
kt = kw*polyval(S1,1)/(dcgain(Wprime)*polyval(T_tilde,1))
T1 = kt*T_tilde;
F = zpk(tf(T1,S1,Ts))

d2 = 0;
rho = 1;

out =  sim('quiz_03_simulink')
figure(2)
plot(out.y.time,out.y.data)
grid on 
xlabel('t')
ylabel('y(t)')
figure(3)
plot(out.u.time,out.u.data)
grid on 
xlabel('t')
ylabel('u(t)')

d2 = 1;
rho = 0;

out =  sim('quiz_03_simulink')
figure(4)
plot(out.y.time,out.y.data)
grid on 
xlabel('t')
ylabel('y(t)')