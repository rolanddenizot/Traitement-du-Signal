clear all;
close all;
clc;

n=1:1:200;
N=200;
fe=1;
Te=1/fe;
t=0:Te:(N-1)*Te;
freq = (-N/2:(N/2-1))*fe/N;

f0=0.25;
r=1;

x1=zeros(1,N);
x1(1)=1;

x2=ones(1,200);

num=[1 0 0];
denum=[1 -2*r*cos(2*pi*f0*Te) r*r];
sys = tf(num, denum, Te);

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

figure(4)
pzplot(sys)
grid on;
zoom on;