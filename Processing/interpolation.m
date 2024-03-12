function [traces]=interpolation(traces,gap)

% Interpolation of missing traces using the mean of neighboring traces
%
% [traces]=interpolation(traces,gap)
%
% Dr. Tina Wunderlich, CAU Kiel 2020/2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces along profile (traces in columns), missing
% traces are nan
% gap: maximum number of traces, which shall be interpolated
%
% Output:
% traces: Matrix with traces using mean of neighboring traces

x=find(~isnan(traces(1,:))); % ~isnan in first sample row = value
trnum=1:length(traces(1,:)); % all trace numbers
miss=find(isnan(traces(1,:))); % number of missing traces

% interpolate everything
ns=length(traces(:,1));
for i=1:ns % for each sample
    traces(i,:)=interp1(x,traces(i,x),trnum,'linear','extrap');
end

% create mask taking into account the gap
mask=NaN(size(trnum));
mask(x)=1; % set given traces to valid
for i=1:length(miss)
    if min(abs(miss(i)-x))<gap/2
        mask(miss(i))=1;
    end
end

% apply mask
mask=repmat(mask,[ns,1]);
traces=traces.*mask;
