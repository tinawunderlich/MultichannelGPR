function [traces]=spectralWhitening(traces,dt,fmin,fmax,alpha)

% Spectral Whitening, based on the description in M. W. Lee (1986)
%
% [traces]=spectralWhitening(traces,dt,fmin,fmax,alpha)
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: Matrix with traces of one channel in columns
% dt: Sampleinterval in ns
% fmin, fmax: minimum and maximum frequency in MHz for bandpass filter
% alpha: Parameter for filtering (if =[], default=0.01)
%
% Output:
% traces: spectrally whitened traces in same order as input 

if isempty(alpha)
    alpha=0.01;
end

% Samplingfrequency
df=1/dt;

% Length of trace
L=length(traces(:,1));

nfft=2^nextpow2(L);

Aw=zeros(nfft,1);
for i=1:length(traces(1,:))
    % autocorrelation:
    xc=xcorr(traces(:,i),traces(:,i),'coeff');
    
    % sum amplitude spectrum of autocorr trace
    Aw=Aw+abs(fft(xc,nfft)./L);
end
Aw=Aw./length(traces(1,:)); % mean amplitude spectrum of autocorrelation trace

% FFT of original traces:
Sw=fft(traces,nfft)./L;

% apply filter:
Swnew=Sw.*Aw.^(alpha-1);

% convert back to time domain:
traces=ifft(Swnew,nfft,'symmetric')./sqrt(alpha); % and scale to get approx. same amplitude level

traces=bandpass_gpr(traces,dt,fmin,fmax); % apply bandpass
