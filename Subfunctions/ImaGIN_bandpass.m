function xf = ImaGIN_bandpass(x, Fs, Fs1, Fp1, Fp2, Fs2)
% Bandpass filter for the signal x.  An acausal fft
% algorithm is applied (i.e. no phase shift). The filter function is
% constructed from a Hamming window (default window used in "fir2" Matlab function).
% to avoid ripples in the frequency response (windowing is a smoothing in frequency domain)
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) frequencies are defined in Hz as
%                  ----------
%                /|         | \
%               / |         |  \
%              /  |         |   \
%             /   |         |    \
%   ----------    |         |     -----------------
%                 |         |
%           Fs1  Fp1       Fp2   Fs2
%
% If NO OUTPUT arguments are assigned the filter function H(f) and
% impulse response are plotted.
%
% OPTIMIZATIONS vs. original:
%   1. Filter H is cached via 'persistent' – recomputed only when parameters
%      or signal length change. Saves fir2+fft on every repeated call.
%   2. Column-wise FFT loop replaced by matrix FFT (fft operates on all
%      columns simultaneously), removing the per-column overhead.
%   3. Unnecessary transpose gymnastics cleaned up.
%
%------------------------------------------------------------------------
% Originally produced by the Helsinki University of Technology,
% Adapted by Mariecito SCHMUCKEN 2001
% (Edited by Dr. Tina Wunderlich, CAU Kiel, 2020,
% tina.wunderlich@ifg.uni-kiel.de)
% (Optimized 2026 – caching + matrix FFT (Claude.ai))
%------------------------------------------------------------------------

% --- Cache filter coefficients across calls ----------------------------
persistent H_cached params_cached

% Ensure column-oriented input
if size(x,1) == 1
    x = x';
end

% Make number of samples EVEN (required by the symmetric FFT trick)
Norig = size(x, 1);
if rem(Norig, 2)
    x = [x; zeros(1, size(x,2))];
end
N  = size(x, 1);
Nh = N / 2;

% Normalised frequency edges
Ns1 = Fs1 / (Fs/2);
Ns2 = Fs2 / (Fs/2);
Np1 = Fp1 / (Fs/2);
Np2 = Fp2 / (Fs/2);

% --- Build / reuse filter -----------------------------------------------
% Key: all parameters that affect H, including signal length N
current_params = [Fs, Fs1, Fp1, Fp2, Fs2, N];
if isempty(params_cached) || ~isequal(current_params, params_cached)
    B            = fir2(N-1, [0 Ns1 Np1 Np2 Ns2 1], [0 0 1 1 0 0]);
    H_cached     = abs(fft(B));   % zero-phase (real, symmetric) filter
    params_cached = current_params;
end
H = H_cached;   % H is a row vector of length N

% --- Apply filter to all columns at once --------------------------------
% fft/ifft work column-wise on matrices; H must broadcast as a column.
xf = real(ifft( fft(x) .* H(:) ));   % H(:) → column → broadcasts over cols
xf = xf(1:Norig, :);

% Trim x back too (needed for the plot branch below)
x = x(1:Norig, :);

% --- Optional diagnostic plots (only when called with no output) --------
if nargout == 0
    IPR = real(ifft(H));
    f   = Fs * (0:Nh-1) / N;
    freqz(IPR, 1, f, Fs);
    figure
    subplot(2,1,1)
    plot(f, H(1:Nh));
    xlim([0 2*Fs2])
    title('Filter function H(f)')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    plot((1:Nh)/Fs, IPR(1:Nh))
    xlim([0 2/Fp1])
    xlabel('Time (sec)')
    ylim([min(IPR) max(IPR)])
    title('Impulse response')
    figure
    subplot(2,1,1), periodogram(x(:,1),   hamming(Norig), 1024, Fs);
    subplot(2,1,2), periodogram(xf(:,1),  hamming(Norig), 1024, Fs);
end

end % ImaGIN_bandpass


% =========================================================================
% LOCAL HELPER: fir2  (no Signal Processing Toolbox required)
% =========================================================================
% Copyright (C) 2000 Paul Kienzle – GPL v2+
% Simplified: ramp handling kept but dead branches removed for clarity.
function b = fir2(n, f, m, grid_n, ramp_n, window)
    if nargin < 3 || nargin > 6
        error('b = fir2(n, f, m [, grid_n [, ramp_n]] [, window])');
    end
    t = length(f);
    if t < 2 || f(1) ~= 0 || f(t) ~= 1 || any(diff(f) < 0)
        error('frequency must be nondecreasing, starting at 0 and ending at 1');
    end
    if t ~= length(m)
        error('frequency and magnitude vectors must be the same length');
    end
    if nargin < 4, grid_n = 512;          end
    if nargin < 5, ramp_n = grid_n / 20;  end

    % Handle window passed in the grid_n / ramp_n slot
    w = [];
    if length(grid_n) > 1, w = grid_n; grid_n = 512;           end
    if length(ramp_n) > 1, w = ramp_n; ramp_n = grid_n / 20;  end
    if nargin < 6,  window = w;              end
    if isempty(window), window = hamming(n+1); end
    if ~isreal(window), window = feval(window, n+1); end
    if length(window) ~= n+1
        error('window must be of length n+1');
    end
    if 2*grid_n < n+1, grid_n = 2^nextpow2(n+1); end

    % Apply ramps to discontinuities
    if ramp_n > 0
        basef = f;  basem = m;
        idx = find(diff(f) == 0);
        f(idx)   = f(idx)   - ramp_n/grid_n/2;
        f(idx+1) = f(idx+1) + ramp_n/grid_n/2;
        idx = find(diff(f) < 0);
        f(idx)   = (basef(idx)   + basef(idx+1)) / 2;
        f(idx+1) = (basef(idx)   + basef(idx+1)) / 2;
        m = interp1(basef, basem, f);
    end

    grid = interp1(f, m, linspace(0, 1, grid_n+1)');
    b    = ifft([grid; grid(grid_n:-1:2)]);
    mid  = (n+1) / 2;
    b    = real([b((2*grid_n - floor(mid) + 1) : 2*grid_n); b(1:ceil(mid))]);

    if size(window, 1) > 1
        b = b .* window;
    else
        b = b' .* window;
    end
    b = b';
end


% =========================================================================
% LOCAL HELPER: hamming  (no Signal Processing Toolbox required)
% =========================================================================
% Copyright (C) 1995-1997 Andreas Weingessel – GPL v2+
function c = hamming(m)
    if nargin ~= 1
        error('hamming(m)');
    end
    if ~(isscalar(m) && m == round(m) && m > 0)
        error('hamming: m must be a positive integer');
    end
    if m == 1
        c = 1;
    else
        c = 0.54 - 0.46 * cos(2 * pi * (0:m-1)' / (m-1));
    end
end