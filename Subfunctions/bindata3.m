function [zout, xgrid, ygrid] = bindata3(z, x, y, xrg, yrg)

% [zout, xgrid, ygrid] = bindata3(z, x, y, xrg, yrg)
%
% based on 2d-code by Patrick Mineault
%    Refs: https://xcorr.net/?p=3326
%          http://www-pord.ucsd.edu/~matlab/bin.htm
% modified/commented by Tina Wunderlich, CAU Kiel 2020-2024
%   tina.wunderlich@ifg.uni-kiel.de
% Claude.ai: Optimized 2026: loop-free accumarray replaces per-sample sparse calls.
%
% x, y     : vectors of length N with x/y coordinates
% z        : traces of size (numsamples x N) corresponding to x/y
% xrg, yrg : vectors with bin edges in x and y direction
% Output:
%   zout   : binned grid (length(yrg)-1) x (length(xrg)-1) x numsamples
%   xgrid, ygrid : coordinate grids (midpoints of bins)
%
% OPTIMIZATIONS vs. original:
%   1. The sample-wise sparse loop is replaced by a single accumarray call
%      that processes all numsamples time slices simultaneously.
%      accumarray(..., @mean) works on the full (bins x numsamples) matrix.
%   2. histc (deprecated) replaced by discretize, which is faster and
%      future-proof.
%   3. Output preallocated as single to match original dtype.

    numsamples = size(z, 1);

    dx = xrg(2) - xrg(1);
    dy = yrg(2) - yrg(1);

    a = xrg + dx/2;
    b = yrg + dy/2;
    [xgrid, ygrid] = meshgrid(a(1:end-1), b(1:end-1));
    [r, c] = size(xgrid);

    % --- Remove coordinates outside the grid ----------------------------
    weg = (x > max(xrg) | x < min(xrg) | y > max(yrg) | y < min(yrg));
    x(weg) = [];
    y(weg) = [];
    z(:, weg) = [];

    % --- Bin assignment -------------------------------------------------
    % discretize returns NaN for out-of-range values; those are already
    % removed above, so no NaN will appear here.
    binsx = discretize(x, xrg);   % column index per trace  (1..c)
    binsy = discretize(y, yrg);   % row    index per trace  (1..r)

    % Clamp last-edge values into last bin (same as original min/max trick)
    binsx = min(max(binsx, 1), c);
    binsy = min(max(binsy, 1), r);

    bins = (binsy - 1) .* c + binsx;   % linear bin index, length N vector

    % --- Loop-free accumulation over ALL time samples at once -----------
    % z is (numsamples x N); we want to average over the N-dimension for
    % each bin independently, producing (numBins x numsamples).
    %
    % accumarray subs must be a column vector; we tile 'bins' for every
    % sample row, and correspondingly expand sample indices.
    %
    % subs(:,1) = bin index  (which spatial bin)
    % subs(:,2) = sample idx (which time sample)

    N    = length(bins);
    bins_col   = bins(:);                          % (N x 1)
    sample_idx = repmat((1:numsamples)', 1, N);    % (numsamples x N)

    subs = [repmat(bins_col, numsamples, 1), ...
            sample_idx(:)];                        % (N*numsamples x 2)

    vals = double(z(:));                           % (N*numsamples x 1), row-major

    % accumarray: output size = (c*r) x numsamples
    zm = accumarray(subs, vals, [c*r, numsamples], @mean, NaN);

    % Reshape to (r x c x numsamples) and cast back to single
    zout = single(permute(reshape(zm, c, r, numsamples), [2 1 3]));

end