function [traces,t]=t0shift(traces,t0,dt)

% [traces,t]=t0correction(traces,t0,dt)
%
% Make t0-correction with constant time shift of t0. Traces
% are padded with zeros to maintain the same length.
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces in columns
% t0: t0 time in same unit as dt
% dt: Sampleinterval in same time unit as t0
%
% Output:
% traces: t0-corrected traces in same order as intraces
% t: new time vector
%

ns=length(traces(:,1));

traces=traces(int16(t0/dt):end,:);    % make t0 correction (cut before t0)
traces(length(traces(:,1))+1:ns,:)=zeros(int16(t0/dt)-1,length(traces(1,:)));   % pad with zeros to original length
t=0:dt:(length(traces(:,1))-1)*dt;