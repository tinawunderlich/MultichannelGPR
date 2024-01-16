function [traces]=applygain(traces,gain)

% [traces]=applygain(traces,gain)
%
% Apply user defined gain curve
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% gain: gainpoints in dB (1. for t=0, last for t=ns*dt) (example: gain=[-20 0 15 25 30])
%
% Output:
% traces: gained traces 
%

ns=length(traces(:,1));

g=interp1(linspace(1,ns,length(gain)),gain,[1:ns]');
traces=traces.*repmat(10.^(g./20),[1 length(traces(1,:))]); % apply gain