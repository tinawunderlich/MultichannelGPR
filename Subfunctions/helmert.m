function dataneu=helmert(data,alt,neu)

%%% dataneu=helmert(data,alt,neu)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: Matrix mit den zu transformierenden Koordinaten in den ersten
% beiden Spalten (x und y) und dann Amplitude in der dritten Spalten (es
% können auch noch mehr Spalten sein, das ist egal)
% alt: alte Koordinaten x,y an mindestens 2 Punkten ([p1x p1y; p2x p2y; ...])
% neu: an den gleichen Punkten die neuen Koordinaten ([p1x p1y; p2x p2y; ...])
%
% dataneu: die gleiche Matrix wie data aber im neuen Koordinatensystem


% Schwerpunkte in beiden Systemen berechnen
xs=mean(alt(:,1));
ys=mean(alt(:,2));
Es=mean(neu(:,1));
Ns=mean(neu(:,2));

% Punktkoordinaten um Schwerpunkte reduzieren
dx=alt(:,1)-xs;
dy=alt(:,2)-ys;
dE=neu(:,1)-Es;
dN=neu(:,2)-Ns;

% Transformationsparameter berechnen
n=length(alt(:,1));
a=(sum(dx.*dE)/n+sum(dy.*dN)/n)/(sum(dx.^2+dy.^2)/n);
b=(sum(dx.*dN)/n-sum(dy.*dE)/n)/(sum(dx.^2+dy.^2)/n);

% Massstab berechnen (sollte ungefähr =1 sein)
q=sqrt(a^2+b^2);

% Punkte aus altem System ins neue transformieren
E0=Es-a*xs+b*ys;
N0=Ns-b*xs-a*ys;

Eneu=E0+a*data(:,1)-b*data(:,2);
Nneu=N0+b*data(:,1)+a*data(:,2);

dataneu=data;
dataneu(:,1)=Eneu;
dataneu(:,2)=Nneu;
