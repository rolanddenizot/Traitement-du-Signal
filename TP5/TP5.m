%% Nettoyage
clc;
close all;
clear all;

%% Sujet
D=0.01;
fe= 96000;
f0 = 5000;
Te=1/fe;

M=D*fe;
t=(1:M)*Te;

x=1+sin(2*pi*f0*t);
y=0*x;
r=0.99;
for i=2:M
    y(i)=r*y(i-1)-(x(i)-x(i-1));
end

figure(1)
plot(1000*t,x,1000*t,y)
legend('|xe(t)|','y(t)')
xlabel('Temps en ms')
grid on;
zoom on;

figure(2)
plot(t,x,t,x-y)
legend('|xe(t)|','|xe(t)-y(t)|')
xlabel('Temps en ms')
grid on;
zoom on;
N=length(t);

figure(4)
Freq=(-N/2:(N/2-1))*fe/N;
plot(Freq, (1/fe)*fftshift(abs(fft(x))), Freq, (1/fe)*fftshift(abs(fft(y))));
xlabel('Fréquence en Hz')
legend('|Xe(f)|','|Y(f)')
grid on;
zoom on;