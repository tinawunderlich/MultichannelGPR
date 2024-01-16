function [traces,flag]=t0corr_thresh(traces,threshold)

% [traces,flag]=t0corr_thresh(traces,threshold)
%
% Make t0-correction via threshold
% and cutting before for all traces
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces in columns
% threshold: amplitude threshold for t0
%
% Output:
% traces: t0-corrected traces in same order as intraces
% flag: vector of length(traces(1,:)) with 0 (no threshold met, original trace) and 1 (threshold met and
% corrected)


flag=zeros(1,length(traces(1,:)));

for c=1:length(traces(1,:)) % loop over all traces
    k=1;
    if threshold<0
        while k<=length(traces(:,c)) && traces(k,c)>threshold
            k=k+1;
        end
    else
        while k<=length(traces(:,c)) && traces(k,c)<threshold
            k=k+1;
        end
    end
    if k>=length(traces(:,c))   % no point found
        traces(:,c)=traces(:,c); % set original trace
        flag(c)=0;
    else
        temp=traces(k:end,c); % cut before threshold is reached
        temp(length(temp)+1:length(temp)+k-1)=zeros(k-1,1);
        traces(:,c)=temp;
        flag(c)=1;
    end
end
