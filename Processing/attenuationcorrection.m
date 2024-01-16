function [traces]=attenuationcorrection(traces,t,sigma,eps)

% function [traces]=attenuationcorrection(traces,t,sigma,eps)
%
% correction of anntenuation via multiplication with exp(alpha*z)
%
% Dr. Tina Wunderlich, CAU Kiel 2020, tina.wunderlich@ifg.uni-kiel.de
%
% Input:
% traces: matrix with traces in columns
% t: time vector in ns
% eps: relative permittivity
% sigma: conductivity in S/m
%
% Output:
% traces: traces in same order with corrected attenuation


% constants
eps0=8.854*1e-12;   % As/Vm
mu0=4*pi*1e-7;  % Vs/Am

% Attenuation coefficient:
alpha=sigma/2*sqrt(mu0/(eps*eps0));

% velocity:
v=1/sqrt(eps0*eps*mu0); % m/s

% depth:
z=(t*1e-9)/2*v; % in m

if ~iscolumn(z)
    z=z';
end

% corrected traces
traces=traces.*repmat(exp(alpha.*z),[1 length(traces(1,:))]);