clear all;
close all;
clc;

load signal2.mat;
s = signal;

f1 = 35000;
f2 = 40000;
fe = 100000;
Te=1/fe;
Ts = 1e-3;

N=length(s);
t=(0:1:N-1)*Ts;

freq=(-N/2:1:(N/2-1))*fe/N;

figure(1)
plot(t,s)
xlabel('Temps en ms')
grid on;
zoom on;
title('s(t)')

figure(2)
plot(freq,abs(fftshift(fft(s))))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('S(f)')

r=0.9;

num1=[1 0 0];
denum1=[1 -2*r*cos(2*pi*f1*Te) r*r];

num2=[1 0 0];
denum2=[1 -2*r*cos(2*pi*f2*Te) r*r];

y1 = filter(num1,denum1,s);
y2 = filter(num2,denum2,s);

figure(3)
plot(freq,abs(fftshift(fft(y1))))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('Y1(f)')

figure(4)
plot(freq,abs(fftshift(fft(y2))))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('Y2(f)')

figure(5)
plot(t,y1)
xlabel('Temps en ms')
grid on;
zoom on;
title('s1(t)')

figure(6)
plot(t,y2)
xlabel('Temps en ms')
grid on;
zoom on;
title('s2(t)')