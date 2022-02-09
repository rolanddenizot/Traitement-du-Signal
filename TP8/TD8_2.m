clear all;
close all;
clc;

N = 200;
fe = 1;
freq = (-N/2:(N/2-1))*fe/N;
Te=1/fe;
t=0:Te:(N-1)*Te;

a = -0.8;
b = -1;

x1=zeros(1,N);
x1(1)=1;

x2=ones(1,200);

num = [1 -b];
denum = [1 -a];

y1 = filter(num,denum,x1);
y2 = filter(num,denum,x2);

figure(1)
plot(freq,abs(fftshift(fft(y1))))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('Réponse fréquentielle')

figure(2)
plot(t,y1)
xlabel('Temps en secondes')
grid on;
zoom on;
title('y1(t)')

figure(3)
plot(t,y2)
xlabel('Temps en s')
grid on;
zoom on;
title('y2(t)')