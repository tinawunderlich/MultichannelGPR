function [traces]=bandpass_gpr(traces,dt,f_start,f_end)

% [traces]=bandpass_gpr(traces,dt,f_start,f_end)
%
% bandpassfiltering of traces, uses function ImaGIN_bandpass of M.
% Schmucken 2001
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% dt: sample interval in ns
% f_start, f_end: frequency interval of bandpass (in MHz)
%
% Output:
% traces: filtered traces 


dts=dt*1e-9;    % sample interval in s

% get corner frequencies
f2=f_start*1e6;
f3=f_end*1e6;
f1=max([1 f2-(f3-f2)/4 f2-(f3-f2)/8]);
f4=(f3+(f3-f2)/4);

ns=length(traces(:,1));

for i=1:length(traces(1,:))
    temp=ImaGIN_bandpass([traces(:,i); zeros(ns,1)],1/dts,f1,f2,f3,f4); % do bandpass and pad trace before
    traces(:,i)=temp(1:ns);
end