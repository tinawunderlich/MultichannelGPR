function [cohe]=coherence(data,t,nwavelength)

% [traces_out]=coherence(traces,t,nwavelength)
%
% Coherence calculation (in&cross) for 3D binned blocks
% based on Trinks & Hinterleitner 2020
%
% Input:
% data: 3d-datacube with traces along third dimension
% t: time vector in ns, same length as data(1,1,:)
% nwavelength: number of wavelengths for correlation calculation
% !!! function assumes 400 MHz antennas and v=10 cm/ns -> wave
% length is v/f=0.25 m
%
% Output:
% cohe: coherence data cube with same size as data
%
% Dr. Tina Wunderlich, CAU Kiel, tina.wunderlich@ifg.uni-kiel.de


dt=t(2)-t(1);
nsamp=round((nwavelength*0.25/0.1)/dt); % number of samples for xcorrelation

mask=~isnan(data(:,:,1));
n=numel(data(:,:,1)); %number of grid cells
siz=size(data(:,:,1));
neigh=[-1 0; 1 0; 0 1; 0 -1]; % neighboring points
rmax=length(data(:,1,1));
cmax=length(data(1,:,1));
cohe=NaN(size(data));
WaitMessage = parfor_wait(n, 'Waitbar', true,'ReportInterval',100);
parfor i=1:n
    
    cohetemp{i}=[];
    
    if mask(i)
        [r,c]=ind2sub(siz,i);
        ind=[r+neigh(:,1) c+neigh(:,2)];
        ind(ind(:,1)<1 | ind(:,1)>rmax | ind(:,2)<1 | ind(:,2)>cmax,:)=[];
        
        % calculate crosscorrelation
        for tt=1:length(t)-nsamp
            corr=zeros(length(ind(:,1)),1);
            for j=1:length(ind(:,1))
                corr(j)=1-empcorrcoeff(permute(data(r,c,tt:tt+nsamp),[3 1 2]),permute(data(ind(j,1),ind(j,2),tt:tt+nsamp),[3 1 2]));
            end
            cohetemp{i}(round(tt+nsamp/2))=sum(corr)/length(corr);
        end
        
    end
    WaitMessage.Send;
end
WaitMessage.Destroy;

for i=1:n
    [r,c]=ind2sub(siz,i);
    if ~isempty(cohetemp{i})
        cohe(r,c,1:length(cohetemp{i}))=cohetemp{i};
    end
end