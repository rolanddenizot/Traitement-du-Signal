% TRAVAIL DE ROLAND DENIZOT ET THOMAS D'ANGLEJAN, ESILV A3, TD-D
%% Nettoyage
clc;
close all;
clear all;

%%% Partie 1: Système Discret "Fun" with sounds
%% 1: Chargement des 3 fichiers
filename1='singing.wav';
filename2='ChurchImpulse.wav';
filename3='DungeonImpulse.wav';
[y1,Fs1]=audioread(filename1);                          % Importation du premier fichier (singing)
[y2,Fs2]=audioread(filename2);                          % Importation du deuxième fichier (ChurchImpulse)
[y3,Fs3]=audioread(filename3);                          % Importation du troisième fichier (DungeonImpulse)

%% 2: Fréquence d'échantillonnage
% Fs1 est la fréquence d'échantillonnage du premier signal
% Fs2 est la fréquence d'échantillonnage du deuxième signal
% Fs3 est la fréquence d'échantillonnage du troisième signal
phrase=['Fréquence déchantillonnage de singing.wav: ',num2str(Fs1),' Hz'];
disp(phrase);
phrase=['Fréquence déchantillonnage de ChurchImpulse.wav: ',num2str(Fs2),' Hz'];
disp(phrase);
phrase=['Fréquence déchantillonnage de DungeonImpulse.wav: ',num2str(Fs3),' Hz'];
disp(phrase);

%% 3: Réponse du système
singingChurch=conv(y1,y2);                                    % Convolution entre le signal d'entrée et le signal Church
singingDungeon=conv(y1,y3);                                   % Convolution entre le signal d'entrée et le signal Dungeon

%% 4: Ecoute des signaux
sound(y1,Fs1);                                                % Ecoute du signal singing (avant convolution)
sound(y2,Fs2);                                                % Ecoute du signal church  (avant convolution)
sound(y3,Fs3);                                                % Ecoute du signal dungeon (avant convolution)
sound(singingChurch,Fs1);                                     % Ecoute de singingChurch  (après convolution)
sound(singingDungeon,Fs1);                                    % Ecoute de singingDungeon (après convolution)

%% 5: Enregitrement des deux résultats sur le PC
audiowrite('singingChurch.wav',singingChurch,Fs1);            % Enregistrement de singingChurch sur le PC
audiowrite('singingDungeon.wav',singingDungeon,Fs1);          % Enregistrement de singingChurch sur le PC

%% 6: Représentation temporelle
% Représentation du premier signal
temps_val=length(singingChurch)/Fs1;
temps=0:1/Fs1:temps_val;                                       % Création du vecteur temps
figure(1);                                                    
plot(temps(1:length(singingChurch)),singingChurch);            % Tracé du signal
xlim([0 5]);                                                   % Bornes de l'axe x
xlabel('Temps en sec')
ylabel('Amplitude du son')
title('Représentation temporelle du signal singingChurch')
grid on;
zoom on;

% Représentation du deuxième signal
temps_val=length(singingDungeon)/Fs1;
temps=0:1/Fs1:temps_val;                                        % Création du vecteur temps
figure(2);
plot(temps(1:length(singingDungeon)),singingDungeon);           % Tracé du signal
xlim([0 5]);                                                    % Bornes de l'axe x
xlabel('Temps en sec')
ylabel('Amplitude du son')
title('Représentation temporelle du signal singingDungeon')
grid on;
zoom on;
 
%% 7: Calcul de la FFT
x4=fft(singingChurch);                  %FFT de singingChurch
x5=fft(singingDungeon);                 %FFT de singingDungeon
 
%% 8: Représentation fréquentielle
% Premier spectre
figure(3);
Xf4=fftshift(abs(x4));                  % Spectre centré
% Vecteur fréquence
N=length(Xf4);                          % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(Fs1/N*I)-Fs1/2;                   % Vecteur de fréquence entre entre -Fs1/2 et +Fs1/2 de taille N
% Représentation
plot(Freq,Xf4);                         % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)| signal singingChurch')
grid on;
zoom on;

% Deuxième spectre
figure(4);
Xf5=fftshift(abs(x5));                  % Spectre centré
% Vecteur fréquence
N=length(Xf5);                          % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(Fs1/N*I)-Fs1/2;                   % Vecteur de fréquence entre entre -Fs1/2 et +Fs1/2 de taille N
% Représentation
plot(Freq,Xf5);                         % Tracé du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)| signal singingDungeon')
grid on;
zoom on;

%%% Partie 2: Quantification et reconditionnement d'un signal audio
%% Question 1
%% 1.a
filename6='extrait.wav';
[y,Fe]=audioread(filename6); % Chargement du fichier
% Fe est la fréquence d'échantillonnage du signal
phrase=['Fréquence déchantillonnage de extrait.wav: ',num2str(Fe),' Hz'];
disp(phrase);
% La fréquence critique au sens de Shannon est: fmax=Fe/2 pouvant être échantillonné sans causer de repliement spectral
fmax=Fe/2;
phrase=['La fréquence critique au sens de Shannon est: ',num2str(fmax),' Hz pouvant être échantillonné sans causer de repliement spectral'];
disp(phrase);

%% 1.b
% Représentation du signal
temps_val=length(y)/Fe;
temps=0:1/Fe:temps_val;                             % Création du vecteur temps
figure(5);
plot(temps(1:length(y)),y)                          % Tracé du signal
xlabel('Temps échantillonné en sec')
ylabel('Amplitude du son')
title('Signal en fonction du temps échantillonné')
grid on;
zoom on;

%% 1.c
bits=audioinfo('extrait.wav');                      % Récupération des infos du signal
disp(bits);
% Le signal est donc quantifié sur 16 Bits
bits=16;

%% 1.d
% On a choisi une fréquence de cet ordre de grandeur car cela correspond à 2
% fois la fréquence maximale qu'un humain peut entendre, il n'y aura donc
% pas de repliement spectral, Shannon est respecté

%% 1.e
nb_echantillons=length(y);
phrase=['Le nombre d échantillons est: ',num2str(nb_echantillons)];
disp(phrase)
temps_total=nb_echantillons/Fe;
phrase=['Le temps total est: ',num2str(temps_total),' s'];
disp(phrase);
% A l'aide de VLC, la valeur est correcte

%% 1.f
taille=nb_echantillons*16/8;                        % Le nombre d'échantillons muliplié par le nombre de bits de quantification et divisé par 8 pour avoir le nombre d'octets
phrase=['La taille est: ',num2str(taille),' o'];
disp(phrase);
taille=taille/1000;
phrase=['La taille est: ',num2str(taille),' Ko'];
disp(phrase);
% Les propriétés du fichier sur l'ordinateur indiquent une taille de
% 1 246 334 octets, ce qui implique une erreur de 0,004%

%% 1.g
sound(y,Fe); % Ecoute du signal

%% Question 2
%% 2.a
precision_theorique=2/(2^bits);         % La précision théorique est l'amplitude max/2^nbbits
phrase=['La précision théorique est d amplitude: ',num2str(precision_theorique)];
disp(phrase);

%% 2.b
precision_pratique=abs(min(y)-max(y))/2^16;
phrase=['La précision pratique est d amplitude: ',num2str(precision_pratique)];
disp(phrase);
% Les valeurs sont du même ordre de grandeur, la valeur pratique est
% cohérente avec la valeur théorique

%% Question 3
%% 3.a
amplitude_min=min(abs(y));           % On determine l'amplitude minimale
amplitude_max=max(y);                % On determine l'amplitude maximale
phrase=['L amplitude minimale est: ',num2str(amplitude_min)];
disp(phrase);
phrase=['L amplitude maximale est: ',num2str(amplitude_max)];
disp(phrase);

%% 3.b
figure(7);
histogram(y,100);               % Création de l'histogramme sur 100 classes
grid on;
zoom on;
% Le signal n'est pas centré en 0, il semble être centré entre 0,05108 et 0,058
% La valeur moyenne du signal à l'aide de l'histogramme semble être 0,05
% environ
val_moyenne=mean(y);            % Calcul de la valeur moyenne du signal
phrase=['La valeur moyenne est: ',num2str(val_moyenne)];
disp(phrase);
% Les valeurs observées à l'oeil nu sur l'histogramme et les valeurs
% calculées sont très proches, donc cohérentes

%% 3.c
ecart_type=std(y);              % Calcul de l'écart-type du signal
phrase=['L ecart-type est: ',num2str(ecart_type)];
disp(phrase);
% L'écart-type est très proche de 0 mais la faible différence explique que
% l'histogramme ne soit pas centré en 0

%% 3.d
% La valeur moyenne étant différente de 0, le signal comporte une
% composante continue, donc le signal n'est pas bien conditionné, nous ne
% tirons pas le meilleur parti du nombre de bits alloués à la
% quantification

%% Question 4
%% 4.a
x=fft(y);                               % Transformation de Fourier rapide de notre signal

%% 4.b
N=length(x);                            % Taille du vecteur temps
I=0:N-1;                                % Même taille que le vecteur temps
Freq=(Fe/N*I)-Fe/2; 

%% 4.c
figure(8);
plot(Freq,x);                           % Représentation fréquentielle du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('Spectre non centré de |X(f)|')
grid on;
zoom on;

%% 4.d
% On constate que le spectre n'est pas symétrique par rapport au point 0,
% il n'est donc pas centré
figure(9);
Xf=fftshift(abs(x)); % Spectre centré
plot(Freq,Xf);       % Représentation fréquentielle du signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('Spectre centré de |X(f)|')
grid on;
zoom on;

%% Question 5
%% 5.a
nouveau_signal=detrend(y); % Création d'un signal centré

%% 5.b
figure(10);
histogram(nouveau_signal,100); % Représentation de l'histogramme sur 100 classes
grid on;
zoom on;
% On remarque que l'histogramme est centré en 0 et la moyenne est égale à
% 0, la composante principale a donc été supprimée. Le signal a bien été
% reconditionné

%% 5.c
% Représentation du signal
temps_val=length(nouveau_signal)/Fe;
temps=0:1/Fe:temps_val;                                 % Création du vecteur temps
figure(11);
plot(temps(1:length(nouveau_signal)),nouveau_signal)    % Représentation du nouveau signal
xlabel('Temps en sec')
ylabel('Amplitude du son')
title('Signal reconditionné en fonction du temps')
grid on;
zoom on;

%% 5.d
Xf_nouveau=fftshift(abs(fft(nouveau_signal)));
N=length(Xf_nouveau);                            % Taille du vecteur temps
I=0:N-1;                                         % Même taille que le vecteur temps
Freq=(Fe/N*I)-Fe/2;
figure(12);
plot(Freq,Xf_nouveau);                           % Représentation fréquentielle du nouveau signal
xlabel('Fréquence en Hz')
ylabel('Représentation spectrale')
title('|X(f)|')
grid on;
zoom on;

%% 5.e
sound(nouveau_signal,Fe);
% La qualité du signal semble meilleure que le signal original, l'audio est
% plus agréable à écouter

%% Question 6
%% 6.a
signal_court=nouveau_signal(1:int32(1.8/temps_total*nb_echantillons)); % Création du signal court qui ne prend que les 1.8 premières secondes du nouveau signal

%% 6.b
sound(signal_court,Fe); % Ecoute du signal

%% 6.c
t_court=0:1/Fe:length(signal_court)/Fe;                 % Création du vecteur temps du signal court

%% 6.d
save('Sauvegarde.mat','signal_court','Fe','t_court');   % Sauvegarde des trois variables demandées dans un fichier matlab

%% 6.e
audiowrite('signal_court.wav',signal_court,Fe);         % Enregistrement de signal_court sur le PC