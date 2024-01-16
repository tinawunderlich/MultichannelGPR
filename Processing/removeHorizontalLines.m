function [traces]=removeHorizontalLines(traces,meanmedian,numtraces)

% [traces]=removeHorizontalLines(traces,meanmedian,numtraces)
%
% removal of horizontal lines by subtracting a mean or median trace
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% optional: meanmedian='mean' or 'median', default is 'mean'
% optional: numtraces: number of traces for moving window over which
% mean/median is calculated. default is 0-> all traces
%
% Output:
% traces: traces in same order with subtracted offset



if nargin==1
    meanmedian='mean';
    numtraces=0;
end
if nargin==2
    numtraces=0;
end
if numtraces==0 % use all traces for mean/median trace
    if strcmp(meanmedian,'mean')
        traces=traces-repmat(mean(traces(:,~isnan(traces(1,:))),2),[1 length(traces(1,:))]);
    else
        traces=traces-repmat(median(traces(:,~isnan(traces(1,:))),2),[1 length(traces(1,:))]);
    end
else % use moving window
    if round((numtraces+1)/2)==numtraces/2+1 % is even...
        numtraces=numtraces+1; % ... make odd
    end
    dtraces=zeros(size(traces));
    if strcmp(meanmedian,'mean')
        dtraces(:,1:(numtraces-1)/2)=traces(:,1:(numtraces-1)/2)-mean(traces(:,~isnan(traces(1,1:(numtraces-1)/2))),2);
        for i=(numtraces-1)/2+1:length(traces(1,:))-(numtraces-1)/2
            dtraces(:,i)=traces(:,i)-mean(traces(:,~isnan(traces(1,i-(numtraces-1)/2:i+(numtraces-1)/2))),2);
        end
        dtraces(:,length(traces(1,:))-(numtraces-1)/2:end)=traces(:,length(traces(1,:))-(numtraces-1)/2:end)-mean(traces(:,~isnan(traces(1,length(traces(1,:))-(numtraces-1)/2:end))),2);
    else
        dtraces(:,1:(numtraces-1)/2)=traces(:,1:(numtraces-1)/2)-median(traces(:,~isnan(traces(1,1:(numtraces-1)/2))),2);
        for i=(numtraces-1)/2+1:length(traces(1,:))-(numtraces-1)/2
            dtraces(:,i)=traces(:,i)-median(traces(:,~isnan(traces(1,i-(numtraces-1)/2:i+(numtraces-1)/2))),2);
        end
        dtraces(:,length(traces(1,:))-(numtraces-1)/2:end)=traces(:,length(traces(1,:))-(numtraces-1)/2:end)-median(traces(:,~isnan(traces(1,length(traces(1,:))-(numtraces-1)/2:end))),2);
    end
    traces=dtraces;
end