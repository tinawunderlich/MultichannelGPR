function [S]=semb(traces,win)

% Subfunction for semblance calculation
%
% [S]=semb(traces,win)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: matrix with traces along columns
% win: window length in time samples (odd!)
%
% Output:
% S: semblance


numtr=length(traces(1,:));
ns=length(traces(:,1));
M=(win-1)/2;

for i=M+1:ns-M
    oben=sum((sum(traces(i-M:i+M,:),2)).^2);
    unten=numtr*sum(sum(traces(i-M:i+M,:).^2,2));
    S(i)=oben/unten;
end
