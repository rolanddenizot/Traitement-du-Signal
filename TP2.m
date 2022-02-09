%% Nettoyage
close all;
clear all;
clc;
%% Base de temps
fe=50;                       % Fréquence d'échantillonnage
Te=1/fe;                     % Période d'échantillonnage
t=-100:Te:100;
%% Base de fréquence
N=length(t);
f=(-N/2:N/2-1)*fe/N;
%% Signal de référence
y=zeros(1,N);
indices=find((t>=-1)&(t<1));
y(indices)=1;                % On met à 1 les points pour lesquels -1<=t<1
%% Calcul du spectre
V=fftshift(abs(fft(y))/fe);  % Module du spectre centré normalisé (calcul numérique)
Y=abs(2*sinc(2*f));          % Module du spectre centré normalisé (calcul analytique) !!!!'sinc' inclut déja pi
%% Tracé
figure;
plot(f,V);
hold all;
plot(f,Y);
legend('V','Y');
grid on
zoom on