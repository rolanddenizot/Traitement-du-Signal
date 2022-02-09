%% Nettoyage
clc;
close all;
clear all;

%% TD 7 Partie C

%% Question 7
load signalMUX.mat
load temps.mat
s=signalMUX;
t=temps;
fe=176400;
Te=1/fe;

%% Question 8
sound(s,fe); % Le signal n'est pas audible car les fréquences sont supérieures à 20kHz

%% Question 9
S=fft(s); % FFT de s(t)
figure(1);
Sf=abs(fftshift(S)); % Spectre centré
% Vecteur fréquence
N=length(Sf); % Taille du vecteur temps
I=0:N-1; % Même taille que le vecteur temps
Freq=(fe/N*I)-fe/2; % Vecteur de fréquence entre entre -Fe/2 et +Fe/2 de taille N
% Représentation
plot(Freq,Sf); % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|S(f)| signalMUX')
grid on;
zoom on;
% Nous constatons clairement que ce signal est en fait consituté de deux
% signaux centrés autour de 20kHz et 35kHz. La bande passante à fc1 est [15000,25000]Hz et à fc2 est [30000,40000]Hz

%% Question 10
%a
fc1=20000;
fc2=35000;
% signal 1
[b1,a1] = butter(3,[fc1-5000,fc1+5000]*2/fe);
sfil1=filter(b1,a1,s); % sfill1 pour "signal filtré 1"
mDemod1 = 2*sfil1;
mDemodf1 = fftshift(fft(mDemod1));
figure(2)
plot(Freq, 1/fe*abs(mDemodf1), Freq,1/fe*abs(mDemodf1));
legend('|S(f)|', 'Spectre du signal demodulé')
xlabel('Fréquence en Hz')
grid on;
zoom on;
% signal 2
[b2,a2] = butter(3,[fc2-5000,fc2+5000]*2/fe);
sfil2=filter(b2,a2,s); % sfill2 pour "signal filtré 2"
mDemod2 =2*sfil2;
mDemodf2 = fftshift(fft(mDemod2));
figure(3)
plot(Freq, 1/fe*abs(mDemodf2), Freq,1/fe*abs(mDemodf2));
legend('|S(f)|', 'Spectre du signal demodulé')
xlabel('Fréquence en Hz')
grid on;
zoom on;

%b
m1=sfil1.*cos(2*pi*fc1*t);
m2=sfil2.*cos(2*pi*fc2*t);

%c
[b3,a3] = butter(3,4000*2/fe);
sfil3=filter(b3,a3,m1); % sfill3 pour "signal filtré 3"
[b4,a4] = butter(3,4000*2/fe);
sfil4=filter(b4,a4,m2); % sfill1 pour "signal filtré 4"

%d
% m1(t) = sfil1(t)/2
% m2(t) = sfil2(t)/2

%e
% Représentation temporelle de m1(t)
N1 = length(sfil3);
t1 = (1:N1)*Te;
figure(4)
plot(t1,sfil3)
title('Signal filtré 3 en fonction du temps')
xlabel('Temps en seconde')
ylabel('Amplitude')
zoom on;
grid on;
% Représentation fréquentielle de m1(t)
Sfil1=abs(fftshift(fft(sfil1)));
figure(5)
plot(Freq, Sfil1);
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
grid on;
zoom on;
title('|Sfil1(f)|')
% Représentation temporelle de m2(t)
N2 = length(sfil4);
t2 = (1:N2)*Te;
figure(6)
plot(t2,sfil4)
title('Signal filtré 4 en fonction du temps')
xlabel('Temps en seconde')
ylabel('Amplitude')
zoom on;
grid on;
% Représentation fréquentielle de m2(t)
Sfil2=abs(fftshift(fft(sfil2)));
figure(7)
plot(Freq, Sfil2);
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
grid on;
zoom on;
title('|Sfil2(f)|')

%f
%%
sound(m1,176400);                        % Ecoute du son
audiowrite('m1(t).wav',m1,176400);       % Sauvegarde du son
%%
sound(m2,176400);                        % Ecoute du son
audiowrite('m2(t).wav',m2,176400);       % Sauvegarde du son

%% Partie 1: Filtre purement récursif
clc;
close all;
clear all;

%% 1)
% Pour que le filtre soit stable, il faut que le dénominateur soit
% différent de 0, donc a!=Z

%% 2)
% y(n)= (1-a)x(n) + ay(n-1)

%% 3)
for a = [-0.8,0.8,0.95,1.2]
    y1 = zeros(1,200); % Crée un vecteur de taille 200 pour initialiser y1(n)
    y2 = zeros(1,200); % Crée un vecteur de taille 200 pour initialiser y2(n)
    y3 = zeros(1,200); % Crée un vecteur de taille 200 pour initialiser y3(n)
    Fe = 1; % Initialisation de Fe
    Te = 1/Fe; % Initialisation de Te
    N = 200; % Initialisation de N
    f = -Fe/2:Fe/N:Fe/2-Fe/N; % Initialisation de f
    z = exp(1j*2*pi*f*Te); % Initialisation de z
    z1 = 1./z; % Initialisation de z^-1

    % Création de l'impulsion de dirac
    x1 = zeros(1,200);
    x1(1) = 1;

    % Calcul de la réponse impulsionnelle y1(n)
    y1(1) = (1-a)*x1(1); % initialisation de la réponse
    for n = 2:200
        y1(n) = (1-a)*x1(n) + a*y1(n-1);
        %disp(y1(n));
    end

    % Création du signal échelon
    x2 = ones(1,200);
    %disp(x2);

    % Calcul de la réponse indicielle y2(n)
    y2(1) = (1-a)*x2(1); % initialisation de la réponse
    for n = 2:200
        y2(n) = (1-a)*x2(n) + a*y2(n-1);
        %disp(y2(n));
    end

    % Module de H(f)
    figure()
    mHF = (1-a)./(1.-z1*a); % Module de H(f)
    N = 200;
	semilogx(f,mHF);
    title('|H(f)| avec échelle logarithmique')
    xlabel('Echelle logarithmique (Hz)')
    ylabel('Amplitude (dB)')
    zoom on;
    grid on;

    % Création du signal sinusoidal
    fs = 0.1;
    x3 = sin((1:N)*2*pi*fs*Te);
    %disp(x3);

    % Calcul de la réponse y3(n) à x3(n)
    y3(1) = (1-a)*x3(1); % initialisation de la réponse
    for n = 2:200
        y3(n) = (1-a)*x3(n) + a*y3(n-1);
        %disp(y3(n));
    end

    % Affichage des sinaux en fonction du temps
    t = (1:N)*Te;
    % Signal y1
    figure()
    plot(t,y1,t,x1)
    title('y1(t) en fonction du temps')
    xlabel('Temps en seconde')
    ylabel('Amplitude')
    zoom on;
    grid on;

    % Signal y2
    figure()
    plot(t,y2,t,x2)
    title('y2(t) en fonction du temps')
    xlabel('Temps en seconde')
    ylabel('Amplitude')
    zoom on;
    grid on;

    % Signal y3
    figure()
    plot(t,y3,t,x3)
    title('y3(t) en fonction du temps')
    xlabel('Temps en seconde')
    ylabel('Amplitude')
    zoom on;
    grid on;

    % Nouvel affichage du module de H(f) en dB avec échelle logarithmique à
    % partir de la réponse y1
    figure()
    semilogx(f, fft(y1))
    title('|H(f)| avec échelle logarithmique à partir de y1')
    xlabel('Fréquence en échelle logarithmique (Hz)')
    zoom on;
    grid on;

    % Poles et des zéros de la fonction de transfert
    numerateur = [1-a 0];
    denominateur = [1 -a];
    pz = tf(numerateur,denominateur,Te);
    figure()
    pzplot(pz)
    title('Poles et zéros de la fonction de transfert')
    zoom on;
    grid on;
end
%% 4)
% Pour a = -0.8: Filtre passe-bas
% Pour a = 0.8: Filtre passe-haut
% Pour a = 0.95: Filtre passe-haut
% Pour a = 1.2: Filtre passe-bas
% Les filtres sont stables pour a = -0.8, a = 0.8 et a = 0.95 car tous les poles
% sont dans le cercle unité, il n'est pas stable pour a = 1.2 car les poles
% sont hors du cercle unité

% Impact de x3(n) sur le passe-haut: Augmente l'amplitude
% Impact de x3(n) sur le passe-bas:

%% Partie 2: Filtre du premier ordre

clc;
close all;
clear all;

for a = [-0.8,0.8]
    for b = [-1,1]
        N = 200; % Initialisation de N
        Fe = 1; % Initialisation de Fe
        Te = 1/Fe; % Initialisation de Te
        f = -Fe/2:Fe/N:Fe/2-Fe/N; % Initialisation de f
        t = 0:Te:(N-1)*Te; % Initialisation de t

        % Initialisation de x1 comme précédemment
        x1 = zeros(1,200);
        x1(1) = 1;

        % Initialisation de x2 comme précédemment
        x2 = ones(1,200);

        numerateur = [1 -b];
        denominateur = [1 -a];

        % Réponse impulsionnelle y1
        y1 = filter(numerateur, denominateur, x1);
        figure()
        plot(t,y1)
        title('y1(t) en fonction du temps')
        xlabel('Temps en seconde')
        ylabel('Amplitude')
        zoom on;
        grid on;

        % Réponse indicielle y2
        y2 = filter(numerateur, denominateur, x2);
        figure()
        plot(t,y2)
        title('y2(t) en fonction du temps')
        xlabel('Temps en seconde')
        ylabel('Amplitude')
        zoom on;
        grid on;

        % Réponse en fréquence
        figure()
        plot(f,abs(fftshift(fft(y1))))
        title('Réponse en fréquence de Y1(f)')
        xlabel('Frequence (Hz)')
        ylabel('Représentation spectrale')
        grid on;
        zoom on;
    end
end

% Pour a = -0.8 et b = -1 : Filtre passe-bande (la bande est beaucoup plus large)
% Pour a = -0.8 et b = 1 : Filtre coupe-bande (la bande est beaucoup plus large)
% Pour a = 0.8 et b = -1 : Filtre passe-bande (la bande est moins large)
% Pour a = 0.8 et b = 1 : Filtre coupe-bande (la bande est moins large)

%% Partie 3: Filtre du second ordre

clc;
close all;
clear all;

N = 200; % Initialisation de N
Fe = 1; % Initialisation de Fe
Te = 1/Fe; % Initialisation de Te
f = -Fe/2:Fe/N:Fe/2-Fe/N; % Initialisation de f
t = 0:Te:(N-1)*Te; % Initialisation de t

for r = [0,0.5,1,1.5,2]
    f0 = 0.1;
    
    % Initialisation de x1 comme précédemment
    x1 = zeros(1,200);
    x1(1) = 1;
    
    % Initialisation de x2 comme précédemment
    x2 = ones(1,200);
    
    numerateur=[1 0 0];
    denominateur=[1 -2*r*cos(2*pi*f0*Te) r*r];
    
    % Réponse impulsionnelle y1
    y1 = filter(numerateur,denominateur,x1);
    figure()
    plot(t,y1)
    title('y1(t) en fonction du temps')
    xlabel('Temps en seconde')
    ylabel('Amplitude')
    grid on;
    zoom on;
    
    % Réponse indicielle y2
    y2 = filter(numerateur, denominateur, x2);
    figure()
    plot(t,y2)
    title('y2(t) en fonction du temps')
    xlabel('Temps en seconde')
    ylabel('Amplitude')
    zoom on;
    grid on;
    
    % Réponse fréquentielle
    figure()
    plot(f,abs(fftshift(fft(y1))))
    title('Réponse en fréquence de Y1(f)')
    xlabel('Frequence (Hz)')
    ylabel('Représentation spectrale')
    grid on;
    zoom on;
    
    % Poles et zéros de la fonction de transfert
    pz = tf(numerateur, denominateur, Te);
    figure()
    pzplot(pz)
    title('Poles et zéros de la fonction de transfert')
    zoom on;
    grid on;
end


%% Partie 4: FSK

%% Avec signal.mat
clc;
clear all;
close all;

%% 1)
load signal.mat
f1 = 40000; % Initialisation de f1
f2 = 20000; % Initialisation de f2
Ts = 0.001; % Initialisation de Ts
Fe = 100000; % Initialisation de Fe
Te = 1/Fe; % Initialisation de Te
s1 = signal;
N = length(s1); % Initialisation de N
f = -Fe/2:Fe/N:Fe/2-Fe/N; % Initialisation de f
t = 0:Te:(N-1)*Te; % Initialisation de f

% Affichage du signal
figure(1)
plot(t,s1)
title('Signal s1 original')
xlabel('Temps en ms')
ylabel('Amplitude')
grid on;
zoom on;

%% 2)
% Représentation fréquentielle
S1 = fft(s1);
figure(2)
plot(f,S1);
title('Représentation fréquentielle de s1')
xlabel('Fréquence (Hz)')
ylabel('Représentation spectrale')
zoom on;
grid on;

%% 3)
r = 0.9;
numerateur1=[1 0 0];
denominateur1=[1 -2*r*cos(2*pi*f1*Te) r*r];
y1 = filter(numerateur1,denominateur1,s1);
figure(3)
plot(f,abs(fftshift(fft(y1))))
title('Réponse en fréquence de Y1(f)')
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
grid on;
zoom on;

numerateur2=[1 0 0];
denominateur2=[1 -2*r*cos(2*pi*f2*Te) r*r];
y2 = filter(numerateur2,denominateur2,s1); % Filtrage du signal
figure(4)
plot(f,abs(fftshift(fft(y2))))
title('Réponse en fréquence de Y2(f)')
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
grid on;
zoom on;

%% 4)
figure(5)
plot(t,y1)
xlabel('Temps en ms')
grid on;
zoom on;
title('Signal s1 filtré')

figure(6)
plot(t,y2)
xlabel('Temps en ms')
grid on;
zoom on;
title('Signal s1 filtré')

%% Avec signal2.mat
clc;
clear all;
close all;

%% 1)
load signal2.mat
f1 = 35000; % Initialisation de f1
f2 = 40000; % Initialisation de f2
Ts = 0.001; % Initialisation de Ts
Fe = 100000; % Initialisation de Fe
Te = 1/Fe; % Initialisation de Te
s2 = signal;
N = length(s2); % Initialisation de N
f = -Fe/2:Fe/N:Fe/2-Fe/N; % Initialisation de f
t = 0:Te:(N-1)*Te; % Initialisation de f

% Affichage du signal
figure(1)
plot(t,s2)
title('Signal s2 original')
xlabel('Temps en ms')
ylabel('Amplitude')
grid on;
zoom on;

%% 2)
% Représentation fréquentielle
S2 = fft(s2);
figure(2)
plot(f,S2);
title('Représentation fréquentielle de s2')
xlabel('Fréquence (Hz)')
ylabel('Représentation spectrale')
zoom on;
grid on;

%% 3)
r = 0.9;
numerateur1=[1 0 0];
denominateur1=[1 -2*r*cos(2*pi*f1*Te) r*r];
y1 = filter(numerateur1,denominateur1,s2); % Filtrage du signal
figure(3)
plot(f,abs(fftshift(fft(y1))))
title('Réponse en fréquence de Y1(f)')
xlabel('Fréquence (Hz)')
ylabel('Représentation spectrale')
grid on;
zoom on;

numerateur2=[1 0 0];
denominateur2=[1 -2*r*cos(2*pi*f2*Te) r*r];
y2 = filter(numerateur2,denominateur2,s2);
figure(4)
plot(f,abs(fftshift(fft(y2))))
title('Réponse en fréquence de Y2(f)')
xlabel('Fréquence (Hz)')
ylabel('Représentation spectrale')
grid on;
zoom on;

%% 4)
figure(5)
plot(t,y1)
xlabel('Temps en ms')
grid on;
zoom on;
title('Signal s2 filtré')

figure(6)
plot(t,y2)
xlabel('Temps en ms')
grid on;
zoom on;
title('Signal s2 filtré')

%% Partie 5: Filtrage d'un bruit aigu dans un signal
clc;
close all;
clear all;

%% Partie 5.1: Etude du signal "ParoleBruit.wav"
[y,Fe]=audioread("Parole_Bruit.wav"); % Fe est la fréquence d'échantillonnage
Te = 1/Fe; % Période d'échantillonnage

% Représentation Temporelle
N = length(y);
t = 0:Te:(N-1)*Te;
figure(1)
plot(t,y)
title('Parole_Bruit en fonction du temps')
xlabel('Temps en seconde')
ylabel('Amplitude')
zoom on;
grid on;
% Représentation Fréquentielle
f = -Fe/2:Fe/N:Fe/2-Fe/N;
Y = abs(fftshift(fft(y)));
figure(2)
semilogy(f,Y)
title('Représentation spectrale de Parole_Bruit')
xlabel('Fréquence (Hz)')
ylabel('Représentation spectrale')
grid on;
zoom on;

% Fréquences indésirables: 800Hz, 4000Hz, 4200Hz, 4500Hz, 5000Hz

%% Partie 5.2: Etude du signal "ParoleBruit.wav"

%% 1)
fc1 = 800; % Première fréquence de coupure
fc2 = 4000; % Deuxième fréquence de coupure
fc3 = 4200; % Troisième fréquence de coupure
fc4 = 4500; % Quatrième fréquence de coupure
fc5 = 5000; % Cinquième fréquence de coupure

%% 2)
[b1,a1] = butter(6,[fc1-25,fc1+25]*2/Fe, 'stop'); % Coupe-bande à la première fréquence
[b2,a2] = butter(6,[fc2-25,fc2+25]*2/Fe, 'stop'); % Coupe-bande à la deuxième fréquence
[b3,a3] = butter(6,[fc3-25,fc3+25]*2/Fe, 'stop'); % Coupe-bande à la troisième fréquence
[b4,a4] = butter(6,[fc4-25,fc4+25]*2/Fe, 'stop'); % Coupe-bande à la quatrième fréquence
[b5,a5] = butter(6,[fc5-25,fc5+25]*2/Fe, 'stop'); % Coupe-bande à la cinquième fréquence

%% 3)
y1 = filter(b1,a1,y); % Premier filtrage
y2 = filter(b2,a2,y1); % Deuxième filtrage
y3 = filter(b3,a3,y2); % Troisième filtrage
y4 = filter(b4,a4,y3); % Quatrième filtrage
yfil = filter(b5,a5,y4); % Cinquième filtrage

%% 4)
yfil = yfil/max(abs(yfil)); % Reconditionnement

%% 5)
% Représentation temporelle
figure(3)
plot(t,y,t,yfil);
title("Représentation temporelle de Parole_Bruit")
xlabel("Temps en seconde")
ylabel("Amplitude")
legend("Avant filtrage", "Après filtrage")
grid on; 
zoom on;
% Représentation fréquentielle
Yfil = abs(fftshift(fft(yfil)));
figure(4)
plot(f,Y,f,Yfil);
title("Représentation fréquentielle de Parole_Bruit")
xlabel("Fréquence (Hz)")
ylabel("Représentation spectrale")
legend("Avant filtrage", "Après filtrage")
grid on; 
zoom on;

%% 6)
sound(yfil,Fe); % Ecoute du son
audiowrite('ParoleFiltre.wav',yfil,Fe) % Sauvegarde du son

%% 7)
fonction_transfert = tf(b1, a1, Te); % Fonction de transfert

%% 8)
figure(5)
bodeplot(fonction_transfert) % Diagramme de Bode du filtre équivalent

%% 9)
figure(6); 
impulse(fonction_transfert/Fe); % Réponse impulsionnelle