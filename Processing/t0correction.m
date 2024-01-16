function [traces,maxind]=t0correction(traces,reference_trace,t0,dt)

% [traces]=t0correction(traces,reference_trace,t0,dt)
%
% Make t0-correction via cross-correlation (alignment to reference_trace)
% and cutting before t0 (picked on reference trace) for all traces. Traces
% are padded with zeros to maintain the same length.
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces in columns
% reference_trace: Reference trace with same length as other traces
% t0: t0 time in same unit as dt, picked on reference trace
% dt: Sampleinterval in same time unit as t0
%
% Output:
% traces: t0-corrected traces in same order as intraces
% maxind: lag of each trace in samples
%
% requires xcorr.m


ns=length(traces(:,1)); % number of samples
maxind=zeros(1,length(traces(1,:)));

for c=1:length(traces(1,:)) % loop over all traces
    % align all traces to reference trace
    [r,lags]=xcorr([zeros(ns,1); reference_trace(1:ns)],[zeros(ns,1); traces(1:ns,c)]);   % Cross correlation
    
    % search for local maximum around 0 lag
    for i=ns*2:ns*4-1
        if r(i)<=r(i+1)
            maxind(c)=lags(i+1);
        else
            break;
        end
    end
    for i=ns*2:-1:2
        if r(i)<=r(i-1)
            maxind(c)=lags(i-1);
        else
            break;
        end
    end
    
    %maxind(c)=lags(find(r(round(length(r)/2)-ns/2:round(length(r)/2)+ns/2)==max(r(round(length(r)/2)-ns/2:round(length(r)/2)+ns/2)))+round(length(r)/2)-ns/2+1); % find maximum correlation index
    trace_long=zeros(2*ns,1);    % make longer dummy trace
    % insert trace at right position
    trace_long(maxind(c)+ns/2+1:maxind(c)+ns/2+ns)=traces(:,c);
    % cut trace to original length, now aligned with reference trace
    traces(:,c)=trace_long(ns/2+1:ns+ns/2);
end

traces=traces(int16(t0/dt):end,:);    % make t0 correction (cut before t0)
traces(length(traces(:,1))+1:length(traces(:,1))+int16(t0/dt)-1,:)=zeros(int16(t0/dt)-1,length(traces(1,:)));   % pad with zeros to original length
