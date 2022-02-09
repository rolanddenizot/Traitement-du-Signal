clear all;
close all;
clc;

load signalMUX.mat;
load temps.mat;

s = signalMUX;
t = temps;

fe = 176400;
Te = 1/fe;

sound(s,fe);

N = length(s);
freq = (-N/2:(N/2-1))*fe/N;
Sf = fftshift(fft(s));

figure(1)
plot(freq, 1/fe*abs(Sf));
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('|S(f)|')

f1=20000;
f2=35000;
[b1, a1] = butter(3,[f1-5000, f1+5000]/(fe/2),'bandpass');
[b2, a2] = butter(3,[f2-5000, f2+5000]/(fe/2), 'bandpass');
s1=filter(b1,a1,s);
s2=filter(b2,a2,s);

m1=s1.*cos(2*pi*f1*t);
m2=s2.*cos(2*pi*f2*t);

[b12, a12] = butter(3,5000/(fe/2));
[b22, a22] = butter(3,f2/(fe/2));
m1=filter(b12,a12,m1);
m2=filter(b22,a22,m2);

N1 = length(m1);
S1 = fftshift(fft(m1));

figure(2)
plot(temps, m1);
xlabel('Temps')
grid on;
zoom on;
title('m1')

figure(3)
plot(freq, 1/fe*abs(S1));
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('|S1(f)|')

N2 = length(m2);
S2 = fftshift(fft(m2));

figure(4)
plot(temps, m2);
xlabel('Temps')
grid on;
zoom on;
title('m2')

figure(5)
plot(freq, 1/fe*abs(S2));
xlabel('Fréquence en Hz')
grid on;
zoom on;
title('|S2(f)|')

sound(m1,fe);
sound(m2,fe);