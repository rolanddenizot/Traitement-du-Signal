clear all;
close all;
clc;

a=1.2;

n=1:1:200;
N=200;
fe=1;
Te=1/fe;
fs=0.1;
t=0:Te:(N-1)*Te;
freq = (-N/2:(N/2-1))*fe/N;
exp=2*i*pi*freq*Te;

num=[1-a 0];
denum=[1 -a];
sys = tf(num, denum, Te);

x1=zeros(1,N);
x1(1)=1;

x2=ones(1,200);

x3 = sin(2*pi*fs*t);

y1=zeros(1,N);
y1(1)=(1-a)*x1(1);
y2=zeros(1,N);
y2(1)=(1-a)*x2(1);
y3=zeros(1,N);
y3(1)=(1-a)*x3(1);

for i = 2:1:200
    y1(i)=(1-a)*x1(i)+a*y1(i-1);
    y2(i)=(1-a)*x2(i)+a*y2(i-1);
    y3(i)=(1-a)*x3(i)+a*y3(i-1);
end

H=(1-a)./(1-a*exp);

figure(1)
semilogx(freq,20*log(abs(H)))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('|H(f)|')

figure(2)
plot(t,y1)
xlabel('Temps en s')
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
plot(t,y3)
xlabel('Temps en s')
grid on;
zoom on;
title('y3(t)')

figure(5)
semilogx(freq,abs(fftshift(fft(y1))))
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('|H(f)|')

figure(6)
pzplot(sys)