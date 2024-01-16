function [traces]=k_highpass(traces,dx,kcutoff)

% [traces]=k_highpass(traces,dx,kcutoff)
%
% horizontal k-filter (highpass)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns (needs more than 24 traces!!)
% dx: spatial sample interval in m
% kcutoff: wavenumber in 1/m
%
% Output:
% traces: filtered traces 




% find nan traces
na=isnan(traces(1,:));
traces(isnan(traces))=0; % set them 0
           

kmax=1/(dx*2);
kmax=kmax-kmax/10000;    % make a little bit smaller than maximum wavenumber
for i=1:length(traces(:,1))
    traces(i,:) = ImaGIN_bandpass(traces(i,:)',1/dx,kcutoff/2,kcutoff,kmax-0.55,kmax)';
end

% set nan-traces back to nan
traces(:,na)=NaN;