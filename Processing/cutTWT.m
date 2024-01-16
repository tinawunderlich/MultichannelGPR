function [traces_out,t,ns]=cutTWT(traces,t,tmax)

%
% Cut TWT of traces
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: time vector in ns
% tmax: time to cut (in ns)
%
% Output:
% traces_out: cut traces
% t: shortened time vector in ns
% ns: new number of samples
%

traces_out=traces(t<=tmax,:);

ns=length(traces_out(:,1));
dt=t(2)-t(1);
t=0:dt:(ns-1)*dt;


