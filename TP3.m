%% Nettoyage
close all;
clear all;
clc;
%% 1
% Variables d'entrée
f1=1000;
f2=2000;
A1=0.5;
A2=0.5; % Attention: |A1| + |A2| <= 1

% a)
fe=44100;           % Fréquence d'échantillonnage
Te=1/fe;            % Période d'échantillonnage
D=2;                % Durée du signal
t=0:Te:D;           % Vecteur du temps entre 0 et 2 sec

% b)
% Signal temporel x(t)
x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t);

% Représentation temporelle
figure(1)
plot(t,x)
xlabel('temps en sec')
ylabel('amplitude du son')
title('x(t)')
grid on;
zoom on;

% c)
% Ecouter le signal, il faut préciser la fréquence d'échantillonnage fe;
sound(x,fe);

% d)
N=length(x);        % Taille du vecteur discret x
X=1/N*fft(x);       % On divise par le nombre de points N pour faire la fft d'un signal périodique

% e)
% Vecteur fréquence
N=length(t);        % Taille du vecteur temps
I=0:(N-1);          % Même taille que le vecteur temps
Freq=(fe/N*I)-fe/2; % Vecteur de fréquence entre -fe/2 et fe/2 de taille N

% f)
% Spectre fréquentiel
Xf=fftshift(abs(x));% Spectre centré

figure(2);
plot(Freq,Xf);
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)|')
grid on;
zoom on;
%% 2
% le fft confirme avec la forme attendue 
% L'oreille humaine entend entre 20 Hz et 20 kHz 
% Quand on dépasse 22.05 kHz, on ne respecte plus Shannon => le spectre n'est plus cohérent avec la théorie.


%% 3
% Quand on dépasse 22.05 kHz, on ne respecte plus Shannon