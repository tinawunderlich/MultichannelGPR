function [traces,shiftsamples]=channelshift(traces,maxshift,shiftsamples)

% [traces,shiftsamples]=channelshift(traces,maxshift,shiftsamples)
%
% Make t0-correction via shifting of complete channels (channel 1 is taken
% as reference channel and all other channels are shifted from -maxshift to
% +maxshift samples and compared to the first channel. The shift with
% maximum correlation is taken for each channel.)
% To use given shifts for each channel (e.g. application of shifts to
% complete area), give shiftsamples as 3rd input.
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Cells for each channel with matrices with traces in columns
% (number of traces per channel can differ!)
% maxshift: Number of samples to shift +/- for testing (optional, default is 40)
% shiftsamples: shift in samples for each channel (optional)
%
% Output:
% traces: t0-corrected traces in same order as intraces
% shiftsamples: shift with maximum correlation for each channel in samples
%
% requires empcorrcoeff.m


ch=length(traces);   % number of channels
for i=1:ch
    if ~isempty(traces{i})
        ntr(i)=length(traces{i}(1,:));    % number of traces per channel
        ns=length(traces{i}(:,1)); % number of samples
    end
end

if nargin==1
    maxshift=40; % Verschiebung in samplen +- (default)
end


if nargin==2  % find best shifts for each channel
    shiftsamples=zeros(ch,1);

    for i=2:ch
        % shift channels
        for j=-maxshift:maxshift    % Verschiebung in samples
            schieb{i}{j+maxshift+1}=mean([zeros(maxshift+j,ntr(i)); traces{i}; zeros(maxshift-j,ntr(i))],2); % shifted mean trace
        end

        % for every shift calculate correlation
        temp=mean([zeros(maxshift,ntr(1)); traces{1}; zeros(maxshift,ntr(1))],2);   % reference: mean 1st channel
        for j=1:length(schieb{2})   
            C(i,j)=empcorrcoeff(temp,schieb{i}{j}); 
        end

        % search for best correation
        bestj(i)=find(max(C(i,:))==C(i,:));
        v=-maxshift:maxshift; % Verschiebung
        shiftsamples(i)=v(bestj(i)); % Anzahl sample von bester verschiebung
    end
end


% given shifts for each channel -> apply 
for i=1:ch
    if ~isempty(traces{i})
        if shiftsamples(i)>0
            traces{i}=[zeros(shiftsamples(i),ntr(i)); traces{i}(1:ns-shiftsamples(i),:)];
        elseif shiftsamples(i)<0
            traces{i}=[traces{i}(-shiftsamples(i)+1:ns,:); zeros(-shiftsamples(i),ntr(i))];
        else % no shift
            traces{i}=traces{i};
        end
    else
        traces{i}=[];
    end
end
    
