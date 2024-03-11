function [traces]=medfilt(traces,m)

% [traces]=medfilt(traces,m)
%
% Apply median filter along individual traces over m samples
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% m: number of samples for moving window
%
% Output:
% traces: filtered traces 
%

n=length(traces(:,1)); % number of samples

% pad data with first/last values
traces=[repmat(traces(1,:),[(m-1)/2 1]); traces; repmat(traces(end,:),[(m-1)/2 1])];

datanew=traces;
for i=1:n-m % go through samples
    datanew(i+(m-1)/2,:)=median(traces(i:i+m-1,:),1); % set median of m values in middle point (for all traces)
end

traces=datanew((m-1)/2+1:(m-1)/2+n,:); % cut padding