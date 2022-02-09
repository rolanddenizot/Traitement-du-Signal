% TRAVAIL DE ROLAND DENIZOT ET THOMAS D'ANGLEJAN, ESILV A3, TD-D
%% Nettoyage
clc;
close all;
clear all;

%% Exercice 2: Tracer les spectres en amplitude des filtres

% Semblable pour les 3 figures
Fe=1;
Te=1/Fe;
N=1000;
f=-Fe/2:Fe/N:Fe/2-Fe/N;
z=exp(1j*2*pi*f*Te);

% Premier filtre
H1=0.5./(1-0.5*z.^-1);
figure(1)
plot(f,20*log10(abs(H1)),f,0*f-3)
legend('Spectre H1 en dB','-3dB')
xlabel('Fréquence f en Hz')
ylabel('|H1(f)| en dB')
grid on;
zoom on;
title('H1 - passe bas');

% Deuxième filtre
H2=0.36./(1+0.64*exp(-4*pi*1j*f*Te));
figure(2)
plot(f,20*log10(abs(H2)),f,0*f-3)
legend('Spectre H2 en dB','-3dB')
xlabel('Fréquence f en Hz')
ylabel('|H2(f)| en dB')
grid on;
zoom on;
title('H2 - passe bande');

% Troisième filtre
H3=0.25-0.5*exp(-1j*2*pi*f*Te)+0.25*exp(-1j*4*pi*f*Te);
figure(3)
plot(f,20*log10(abs(H3)),f,0*f-3)
legend('Spectre H3 en dB','-3dB')
xlabel('Fréquence f en Hz')
ylabel('|H3(f)| en dB')
grid on;
zoom on;
title('H3 - passe bande');

%% TP1: Filtre Numérique et Matlab
H1=tf(0.5,[1 -0.5]);
figure(3);
bodeplot(H1);                       % Diagramme de bode du premier filtre
H2=tf(0.36,[1 0.64]);
figure(4);
bodeplot(H2);                       % Diagramme de bode du deuxième filtre
H3=tf([0.25 -0.5 0.25]);
figure(5);
bodeplot(H3);                       % Diagramme de bode du troisième filtre
figure(10);
impulse(H1);            % Réponse impulsionnelle du premier filtre
figure(11);
impulse(H2);            % Réponse impulsionnelle du deuxième filtre
figure(12);
impulse(H3);            % Réponse impulsionnelle du troisième filtre

%% TP2: Filtrage d'un bruit aigu dans un signal
%% A - Etude du signal "Mozart_Bruit.wav"

% Question 1
[y1,Fs1]=audioread("Mozart_Bruit.wav"); 
disp(Fs1) % Fs1 est la fréquence d'échantillonnage du signal
periode=1/Fs1; % La période est 1/Fs1

% Question 2
% Représentation temporelle
temps_val=length(y1)/Fs1;
temps=0:1/Fs1:temps_val;                                       % Création du vecteur temps
figure(6);                                                    
plot(temps(1:length(y1)),y1);                                  % Tracé du signal
xlabel('Temps en sec')
ylabel('Amplitude du son')
title('Représentation temporelle du signal')
grid on;
zoom on;

% Représentation fréquentielle
x=fft(y1);                              % FFT de mozart
figure(7);
Xf=fftshift(abs(x));                    % Spectre centré
% Vecteur fréquence
N=length(Xf);                           % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(Fs1/N*I)-Fs1/2;                   % Vecteur de fréquence entre entre -Fs1/2 et +Fs1/2 de taille N
% Représentation
plot(Freq,Xf);                          % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)| signal mozart_bruit')
grid on;
zoom on;

% La fréquence indésirable est entre 5650 et 5700 Hz

%% B - Filtrage numérique en post-traitement (filtre idéal non-causal)
% Il faut un filtre passe-bas afin de supprimer les fréquences hautes,
% correspondant aux aigus
H=0;
son_filtre=conv(y1,H);          % Génération du son
sound(son_filtre,Fs1);          % Ecoute du son
audiowrite('Mozart_Filtre.wav',son_filtre,Fs1); % Sauvegarde du son

%% C - Filtrage Numérique en temps réel (filtre IIR causal)
% Ce filtre devra couper la bande [5650 5700] Hz
[b,a]=butter(1,5650/(5700/2));          % Coefficients du filtre
signalFiltre=filter(b,a,y1);            % Filtrage du signal d'origine

% Représentation temporelle
temps_val=length(signalFiltre)/Fs1;
temps=0:1/Fs1:temps_val;                                       % Création du vecteur temps
figure(8);                                                    
plot(temps(1:length(signalFiltre)),signalFiltre);              % Tracé du signal
xlabel('Temps en sec')
ylabel('Amplitude du son')
title('Représentation temporelle du signal MozartFiltre')
grid on;
zoom on;

% Représentation fréquentielle
x=fft(signalFiltre);                    % FFT de mozartFiltre
figure(9);
Xf=fftshift(abs(x));                    % Spectre centré
% Vecteur fréquence
N=length(Xf);                           % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(Fs1/N*I)-Fs1/2;                   % Vecteur de fréquence entre entre -Fs1/2 et +Fs1/2 de taille N
% Représentation
plot(Freq,Xf);                          % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)| signal MozartFiltre')
grid on;
zoom on;

sound(signalFiltre,Fs1);                % Ecoute du son
audiowrite('MozartFiltre.wav',son_filtre,Fs1); % Sauvegarde du son

fonction_transfert=tf([b,a]);           % Fonction de transfert
figure(13);
bodeplot(fonction_transfert);           % Diagramme de Bode du filtre équivalent
figure(14);
impulse(fonction_transfert);            % Réponse impulsionnelle du filtre équivalent
