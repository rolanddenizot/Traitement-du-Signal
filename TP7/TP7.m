%% Nettoyage
clc;
close all;
clear all;

%% Partie C
%% Question 7
load signalMUX.mat
load temps.mat
s=signalMUX;
t=temps;
fe=176400;
Te=1/fe;

%% Question 8
sound(s,fe);                            % Le signal n'est pas audible car il contient plusieurs signaux, c'est une porteuse.

%% Question 9
S=fft(s);                               % FFT de s(t)
figure(1);
Sf=fftshift(abs(S));                    % Spectre centré
% Vecteur fréquence
N=length(Sf);                           % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(fe/N*I)-fe/2;                     % Vecteur de fréquence entre entre -Fe/2 et +Fe/2 de taille N
% Représentation
plot(Freq,Sf);                          % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|S(f)| signalMUX')
grid on;
zoom on;
% Nous constatons clairement que ce signal est en fait consituté de deux
% signaux centrés autour de 20kHz et 35kHz. La bande passante à fc1 est [15000,25000]Hz et à fc2 est [30000,40000]Hz

%% Question 10
%a
[b,a] = butter(3,5000*2/fe);
sfil=filter(b,a,s);
mDemod =2*sfil;
mDemodf = fftshift(fft(mDemod));
figure(2)
plot(Freq, 1/fe*abs(mDemodf), Freq,1/fe*abs(mDemodf));
legend('|S(f)|', 'Spectre du signal demodulé')
xlabel('Fréquence en Hz')
grid on;
zoom on;

%b
fc1=20000;
fc2=35000;
s1=s*cos(2*pi*fc1*t);
s2=s*cos(2*pi*fc2*t);

%c

%d

%e

%f
sound(mi,176400);                        % Ecoute du son
audiowrite('mi(t).wav',mi,176400);       % Sauvegarde du son