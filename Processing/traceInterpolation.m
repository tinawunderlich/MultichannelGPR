function [traces]=traceInterpolation(traces,minfactor,maxfactor)

% Replacement of anomalous traces: Strange traces are removed and
% interpolated with the mean of neighboring traces (for single bad traces)
% or with linear interpolation of neighboring traces (for several bad traces)
%
% [traces]=traceInterpolation(traces,minfactor,maxfactor)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces along profile (traces in columns)
% minfactor: remove trace if mean(abs(trace)) is smaller than 1/minfactor*(mean of abs(all traces))
% maxfactor: remove trace if mean(abs(trace)) is higher than maxfactor*(mean of abs(all traces))
%
% Output:
% traces: Matrix with traces, where anomalous traces have been replaced by
% interpolation of neighboring traces


if nargin==1
    % use default values
    minfactor=3;
    maxfactor=3;
end


%mittel=mean(abs(traces(:))); % mean of abs(all traces)

mitteltrace=mean(abs(traces),1); % mean of abs(traces) for each trace

[~,ispikes] = hampel(mitteltrace);
%fprintf("- Number of bad traces found and removed: %04d\n",sum(ispikes));


%chunks=findchunks(mitteltrace>=mittel*maxfactor | mitteltrace<=mittel/minfactor);
chunks=findchunks(ispikes);

if ~isempty(chunks)
    for i=1:length(chunks(:,1))
        if chunks(i,1)==chunks(i,2) % only one trace
            if chunks(i,1)>1 && chunks(i,1)<length(traces(1,:)) % in the middle
                traces(:,chunks(i,1))=mean(traces(:,[chunks(i,1)-1 chunks(i,1)+1]),2); % mean of neighboring traces
            elseif chunks(i,1)==1 % at beginning
                traces(:,chunks(i,1))=traces(:,chunks(i,1)+1);
            elseif chunks(i,1)==length(traces(1,:)) % at end
                traces(:,chunks(i,1))=traces(:,chunks(i,1)-1);
            end
        else
            if chunks(i,1)>1 && chunks(i,2)<length(traces(1,:)) % in the middle
                for j=1:length(traces(:,1))
                    traces(j,chunks(i,1):chunks(i,2))=interp1([chunks(i,1)-1 chunks(i,2)+1],traces(j,[chunks(i,1)-1 chunks(i,2)+1]),chunks(i,1):chunks(i,2));
                end
            elseif chunks(i,1)==1 % at beginning
                traces(:,chunks(i,1):chunks(i,2))=repmat(traces(:,chunks(i,2)+1),[1 chunks(i,2)-chunks(i,1)+1]);
            elseif chunks(i,2)==length(traces(1,:)) % at end
                traces(:,chunks(i,1):chunks(i,2))=repmat(traces(:,chunks(i,1)-1),[1 chunks(i,2)-chunks(i,1)+1]);
            end
        end
    end
end

end

function chunks=findchunks(ind)
% find chunks in ind (=blocks of following traces that have to be replaced)
chunks=[];
flag=0;
for ii=1:length(ind)-1
    if ii==2 && ind(ii)==0 && ind(ii-1)==1
        flag=flag+1;
        chunks(flag-1,2)=ii-1; % End of interval, if only first trace is bad
    end
    if (ind(ii)==0 && ind(ii+1)>0)
        flag=flag+1;
        chunks(flag,1)=ii+1;  % start of interval with data
        if ii==length(ind)-1
            chunks(flag,2)=ii+1; % set end of interval
        end
    elseif (ii==1 && ind(ii)>0) % start of interval for first trace
        flag=flag+1;
        chunks(flag,1)=ii;  % start of interval with data
    elseif (ind(ii)==1 && ind(ii+1)==0)
        flag=flag+1;
        chunks(flag-1,2)=ii;    % end of interval with data
    elseif (ii+1==length(ind) && ind(ii+1)>0) % end of line
        flag=flag+1;
        chunks(flag-1,2)=ii+1;    % end of interval with data
    end
end
if ~isempty(chunks)
    chunks(chunks(:,1)==0 & chunks(:,2)==0,:)=[];
end
end