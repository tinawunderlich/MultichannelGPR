function [data]=envelope(data)


%%% Envelope of each trace in 3D-block
%
% [data]=envelope(data)
%
% Dr. Tina Wunderlich, CAU Kiel 2023, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% data: 3D-block
%
% Output:
% data: envelope of each trace, same size as data


for i=1:length(data(:,1,1))
    for j=1:length(data(1,:,1))
        if ~(all(isnan(data(i,j,:))))
            % take only part with data (if topography corrected/migrated
            % data)
            trace=permute(data(i,j,:),[3 1 2]);
            ind=find(~isnan(trace));
            h=hilbert(trace(ind));

            data(i,j,ind)=sqrt(trace(ind).^2+permute(imag(h)',[2 3 1]).^2);
        end
    end
end

