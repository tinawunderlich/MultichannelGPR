function [traces]=interpolation(traces,gap)

% Interpolation of missing traces using the mean of neighboring traces
%
% [traces]=interpolation(traces,gap)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces along profile (traces in columns), missing
% traces are nan
% gap: maximum number of traces, which shall be interpolated
%
% Output:
% traces: Matrix with traces using mean of neighboring traces




for i=2:length(traces(1,:)) % for all traces
    a=1;
    while a<=gap && i+a<length(traces(1,:))
        if ~isnan(traces(1,i-1)) && isnan(traces(1,i)) && ~isnan(traces(1,i+a))
            traces(:,i:i+a-1)=repmat(mean(traces(:,[i-1 i+a]),2),[1 a]);
            break;
        end
        a=a+1;
    end
end