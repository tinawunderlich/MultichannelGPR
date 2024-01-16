function [f,absAmps]=makeAmpspec(t,traces)

%
% Calculate amplitude spectrum from traces
%
% Dr. Tina Wunderlich, CAU Kiel 2021, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% traces: Matrix with traces in columns
% t: time vector in ns
%
% Output:
% f: frequency vector in MHz
% absAmps: absolute amplitude spectra of traces in columns
%

% Samplinginterval
dt=t(2)-t(1);

% Samplingfrequency
df=1/dt;

% lenght of time series
L=length(t);

% elongated length
nfft=2^nextpow2(L);

% FFT of traces padded with zeros
Y=fft(traces,nfft)/L;

% frequency vector
f=df/2*linspace(0,1,nfft/2+1).*1e3; % in MHz

% Amplitude spektra
absAmps=2*abs(Y(1:nfft/2+1,:));
