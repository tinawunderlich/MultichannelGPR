function [out]=schnittpunkt(koord1,koord2)

%%%% 2 Geraden werden durch jeweils 2 Punkte beschrieben, Programm rechnet
%%%% Schnittpunkt aus
%
% [out]=schnittpunkt(koord1,koord2)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% koord1=[Rw1 Hw1; Rw2 Hw2]
% koord2=[...]
%
% Output:
% out=[Rw_neu Hw_neu]



% [out]=schnittpunkt(koord1,koord2)
% koord1=[Rw1 Hw1; Rw2 Hw2], koord2=[...], out=[Rw_neu Hw_neu]

% % Zur Stabilisierung zentrieren
% mittel_Rw=mean([koord1(:,1); koord2(:,1)]);
% mittel_Hw=mean([koord1(:,2); koord2(:,2)]);
% 
% 
% % Gerade 1 - Gleichung bestimmen
% p1=polyfit([koord1(:,1)]-mittel_Rw,[koord1(:,2)]-mittel_Hw,1);
% 
% % Gerade 2 - Gleichung bestimmen
% p2=polyfit([koord2(:,1)]-mittel_Rw,[koord2(:,2)]-mittel_Hw,1);
% 
% 
% % neuer x-Wert
% out(1)=(p2(2)+p1(2))/(p1(1)+p2(1))+mittel_Rw;
% 
% % neuer y-Wert
% out(2)=polyval(p1,out(1)-mittel_Rw)+mittel_Hw;


flag=0;
% falls eine senkrechte Gerade -> umdrehen
if (koord1(1,1)==koord1(2,1) && koord2(1,2)~=koord2(2,2)) || (koord2(1,1)==koord2(2,1) && koord1(1,2)~=koord1(2,2))
    temp=koord1(:,1);
    koord1(:,1)=koord1(:,2);
    koord1(:,2)=temp;
    temp=koord2(:,1);
    koord2(:,1)=koord2(:,2);
    koord2(:,2)=temp;
    flag=1;
end

% falls beide Geraden senkrecht zu einander und genau in x- und y-Richtung
if koord1(1,1)==koord1(2,1) && koord2(1,2)==koord2(2,2)
    out=[koord1(1,1) koord2(1,2)];
    flag=2;
elseif koord2(1,1)==koord2(2,1) && koord1(1,2)==koord1(2,2)
    out=[koord2(1,1) koord1(1,2)];
    flag=2;
end
    
if flag~=2
    % Steigung Gerade 1
    a1=(koord1(2,2)-koord1(1,2))/(koord1(2,1)-koord1(1,1));
    %Achsenabschnitt Gerade 1
    b1=-a1*koord1(1,1)+koord1(1,2);

    % Steigung Gerade 2
    a2=(koord2(2,2)-koord2(1,2))/(koord2(2,1)-koord2(1,1));
    %Achsenabschnitt Gerade 2
    b2=-a2*koord2(1,1)+koord2(1,2);

    % x-Wert Schnittpunkt
    xs=(b2-b1)/(a1-a2);

    % y-Wert Schnittpunkt
    ys=a1*xs+b1;

    if flag==1
        out=[ys xs];
    else
        % Output
        out=[xs ys];
    end
end