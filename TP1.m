%% Nettoyage
close all;
clear all;
clc;
%% Paramètres de la simulation
T0=1;               % Période du signal
alpha=0.3;          % Rapport cyclique
f0=1/T0;            % Fréquence du signal
w0=2*pi*f0;         % Pulsation du signal
M=10;               % Nombre d'harmoniques à utiliser
%% Base de temps
fe=10e3;
Te=1/fe;
t=-2:Te:2-Te;
%% Signal de référence
x=zeros(1,length(t)/4);
x(1:alpha*length(t)/4)=1;
x=[x,x,x,x];        % Le signal sur l'intervalle désiré est composé de répétitions du motif de base
%% Décomposition de Fourier
xm=zeros(1,length(t));
xm=xm+alpha;
for n=1:M
    xm=xm+(sin(2*n*alpha*pi)*cos(n*w0*t)+(1-cos(2*n*alpha*pi))*sin(n*w0*t))/(n*pi);
end;
%% Puissance sur le spectre
p=alpha^2;
for n=1:M
    p=p+0.5*(((sin(2*n*alpha*pi))/(n*pi))^2+(((1-cos(2*n*alpha*pi)))/(n*pi))^2);
end;
disp(p);
%% Tracé
figure;
plot(t,x);
hold all;
plot(t,xm);