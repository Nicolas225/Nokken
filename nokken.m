% Gegevens :	De te ontwerpen nok moet volgende heffing kunnen realiseren :
% 
% van 10� tot 100� : +35 mm
% van 180� tot 210� : -20 mm
% van 210� tot 280� : -15 mm
% De equivalente massa en dempingsconstante van de volger (en bijhorende onderdelen) worden respectievelijk geschat op 0.40 kg en 0.07, terwijl het mechanisme volgende statische krachten moet uitoefenen:
% 
% van 10� tot 190� : een konstante trek-kracht van 230 N
% van 190� tot 210� : een konstante druk-kracht van 550 N
% van 210� tot 250� : een lineair afnemende trek-kracht van 300 N tot 0 N
% De gevraagde cyclustijd voor de bewerking uitgevoerd door de volger is 0.05 sec.

%% Hefwet

% Berekening van de tweede cyclo�de (= segmenten 2 en 3)


L = 35;
t = 30;


C3 = @(b) (15 - L*(1-(t/b)+(1/pi)*sin(pi*(t/b))));
C6 = @(b) (15 - L*(1-(t/b)+(1/(2*pi))*sin(2*pi*(t/b))));

optim_options = optimset('Display','off','Tolfun',1e-13);

[b3, fval3, exitflag3]=fsolve(C3,[45]',optim_options);
[b6, fval6, exitflag6]=fsolve(C6,[60]',optim_options);
disp(b3 + 180);
disp(b6 + 180);

%% 

out = load('params.mat');

