clear all;
c       = physconst('LightSpeed');

%% --- radar timing params (edit freely; nt2s fitting is generic) ---
N_chirp = 512;
fc0     = 75.9474e9;
tr0     = (60*54) * (1/54e6);   % example B: ~60us
deltaB  = 3e9/512;
% Example A (older):
% tr0     = (48*54) * (1/54e6);
% deltaB  = 1.364501953e6;

v       = 10;
R       = 0;

nX0     = db2mag(20);
nAmp    = db2mag(-57);

tScale  = 1/54e6;
% tScale  = 1/756e6;

%% Fitting controls
Nseg_max   = 10;   % allow up to this many whole-period segments
force_nseg = 0;    % 0: auto (95% SSE gain elbow); >0: exact segment count
fit_stride = 1;    % candidate boundary stride (1=exact; 2/4=faster approx)
% Optional manual lock (skips DP):
use_manual_segs = 0;
manual_segs = [1, 2; 280, 3];   % [start_chirp_1based, period]

%% Precompute ideal nt3 / quantized nt2 (independent of velocity)
for n = 1:N_chirp
    t3(n)  = fc0*tr0/(fc0+(n-1)*deltaB);
    t2(n)  = round(t3(n)/tScale)*tScale;
    t2(n)  = tr0;
    st3(n) = n * fc0*tr0/(fc0+(n-1)*deltaB);
    if n > 1
        nt3(n) = st3(n) - st3(n-1);
    else
        nt3(n) = tr0;
    end
    nt2(n) = round(nt3(n)/tScale)*tScale;
end

%% Detect plateau lengths of ideal nt2
dx = diff(nt2);
jump_indices = find(abs(dx) > tScale/2);
plateau_starts = [1, jump_indices + 1];
plateau_ends   = [jump_indices, N_chirp];
detected_lengths = plateau_ends - plateau_starts + 1;
plateau_values   = nt2(plateau_starts);

fprintf('\n=== nt2 plateau lengths (ideal quantized) ===\n');
fprintf('tScale=%.6f ns, plateaus=%d, lengths unique=%s\n', ...
    tScale*1e9, numel(detected_lengths), mat2str(unique(detected_lengths)));
fprintf('length sequence: %s\n', mat2str(detected_lengths));

inner_lengths = detected_lengths;
if numel(inner_lengths) >= 3
    inner_lengths = detected_lengths(2:end-1);
end
approx_period = round(median(inner_lengths));
period_theory = tScale / mean(abs(diff(nt3)));
fprintf('median plateau length=%d, theory period=%.3f\n', approx_period, period_theory);

%% Build nt2s with a universal piecewise-period fit
% Method: DP minimizes SSE vs ideal nt2 levels, using up to Nseg_max
% segments. Each segment has one integer period P (drop 1 level every P chirps).
% Period candidates come from local theory P(n)=tScale/|d nt3|, so the
% same code works when params change (e.g. periods near 11 or near 2).
levels = round(nt2 / tScale);
L0 = levels(1);
d_nt3 = abs(diff(nt3));
d_nt3(d_nt3 == 0) = tScale;   % avoid /0
P_th = tScale ./ d_nt3;
Pmin = max(1, floor(min(P_th)));
Pmax = max(Pmin, ceil(max(P_th)));
P_list = Pmin:Pmax;

fprintf('\n=== Universal nt2s fit (DP on level SSE) ===\n');
fprintf('theory P range [%.3f, %.3f] -> candidate periods %s\n', ...
    min(P_th), max(P_th), mat2str(P_list));

if use_manual_segs
    seg_table = manual_segs;
    best_lv = build_level_segs(N_chirp, L0, seg_table);
    nseg_used = size(seg_table, 1);
else
    fprintf('DP fitting (N=%d, Nseg_max=%d, stride=%d) ...\n', N_chirp, Nseg_max, fit_stride);
    tic;
    [seg_table, best_lv, nseg_used, seg_costs] = fit_nt2s_dp( ...
        levels, P_list, Nseg_max, force_nseg, fit_stride);
    fprintf('done in %.2fs\n', toc);
    fprintf('SSE by segment count: %s\n', mat2str(seg_costs, 4));
end

nt2s = best_lv * tScale;

fprintf('used %d segment(s) [start_chirp, period]:\n', nseg_used);
fprintf('  (each segment resets the chirp counter; adjacent same-P is OK)\n');
fprintf('  %6s %8s\n', 'start', 'period');
for iseg = 1:size(seg_table, 1)
    fprintf('  %6d %8d\n', seg_table(iseg, 1), seg_table(iseg, 2));
end
fprintf('level RMSE=%.4f, max|err|=%d, end %d (ideal %d)\n', ...
    sqrt(mean((double(levels) - double(best_lv)).^2)), ...
    max(abs(levels - best_lv)), best_lv(end), levels(end));
fprintf('cumsum err = %.6f us\n', (sum(nt2s) - sum(nt2)) * 1e6);

dxs = diff(nt2s);
jump_s = find(abs(dxs) > tScale/2);
len_s = diff([0, jump_s, N_chirp]);
fprintf('nt2s plateau lengths: %s\n', mat2str(len_s));

figure(2);
subplot(311);
stem(detected_lengths, 'filled');
xlabel('plateau index'); ylabel('length');
title(sprintf('nt2 ideal lengths (median=%d, theory=%.2f)', approx_period, period_theory));
grid on;

subplot(312);
stem(len_s, 'filled');
xlabel('plateau index'); ylabel('length');
title(sprintf('nt2s regularized (%d whole-period segments)', nseg_used));
grid on;

subplot(313);
stairs(1:N_chirp, nt2*1e6, 'b'); hold on;
stairs(1:N_chirp, nt2s*1e6, 'r');
yl = [min([nt2, nt2s]), max([nt2, nt2s])] * 1e6;
for iseg = 2:size(seg_table, 1)
    plot([seg_table(iseg, 1), seg_table(iseg, 1)], yl, 'k--');
end
hold off;
xlabel('chirp index n'); ylabel('us');
legend('nt2', 'nt2s', 'Location', 'best');
title('nt2 vs nt2s');
grid on;

%% Velocity sweep: compare nt2 / nt2s / nt3
for v = 0:10:100
    for n = 1:N_chirp
        x0(n)   = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(t2(1:n))) / c);
        x1(n)   = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(t3(1:n))) / c);
        nx0(n)  = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(nt2(1:n))) / c);
        nx0s(n) = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(nt2s(1:n))) / c);
        nx1(n)  = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(nt3(1:n))) / c);
    end
    ny0  = nx0  + nAmp*randn(size(nx0));
    ny0s = nx0s + nAmp*randn(size(nx0s));
    ny1  = nx1  + nAmp*randn(size(nx1));

    figure(1);
    subplot(121);
    plot(nt2); hold on; plot(nt2s); plot(nt3); hold off;
    legend('nt2', 'nt2s', 'nt3', 'Location', 'best');
    title(sprintf('Tr0=%.1fus, N=%d, segs=%d', tr0*1e6, N_chirp, nseg_used));

    subplot(122);
    w = chebwin(N_chirp).';
    plot(mag2db(abs(fft(ny0 .*w, N_chirp)))); hold on;
    plot(mag2db(abs(fft(ny0s.*w, N_chirp))), '-.');
    plot(mag2db(abs(fft(ny1 .*w, N_chirp))), '-o');
    hold off;
    legend('nt2', 'nt2s', 'nt3', 'Location', 'best');
    grid on; grid minor; ylim([-40 80]);
    title(sprintf('scale=%dMHz, vel=%dm/s', 1/tScale/1e6, v));
    pause(0.3);
end

%% ===== helpers =====

function lv = build_level_segs(N, L0, segs)
% Build staircase levels from segment table [start_1based, period].
    lv = zeros(1, N);
    cur = L0;
    nSeg = size(segs, 1);
    for i = 1:nSeg
        s = segs(i, 1);
        P = segs(i, 2);
        if i < nSeg
            e = segs(i+1, 1) - 1;
        else
            e = N;
        end
        seg_base = cur;
        for n = s:e
            lv(n) = seg_base - floor((n - s) / P);
        end
        cur = lv(e);
    end
end

function [seg_table, best_lv, nseg_used, costs] = fit_nt2s_dp(levels, P_list, Nseg_max, force_nseg, fit_stride)
% DP: min SSE(ideal levels, reconstructed levels) with <= Nseg_max segments.
% Each segment uses one integer period from P_list.
    if nargin < 5 || isempty(fit_stride)
        fit_stride = 1;
    end
    fit_stride = max(1, round(fit_stride));

    N = numel(levels);
    L0 = levels(1);
    K = max(1, Nseg_max);
    nP = numel(P_list);
    INF = 1e100;

    % dp(k+1, i+1): best SSE covering chirps 1..i with k segments
    dp = INF * ones(K+1, N+1);
    endL = zeros(K+1, N+1);
    prev = -ones(K+1, N+1);
    prevP = zeros(K+1, N+1);
    dp(1, 1) = 0;
    endL(1, 1) = L0;

    for k = 1:K
        kr = k + 1;
        kprev = k;
        for i = 1:N
            ic = i + 1;
            j_list = 0:fit_stride:i-1;
            if j_list(end) ~= i-1
                j_list = [j_list, i-1]; %#ok<AGROW>
            end
            for j = j_list
                jc = j + 1;
                if dp(kprev, jc) >= INF
                    continue;
                end
                base = endL(kprev, jc);
                ideal = double(levels(j+1:i));
                Lseg = i - j;
                idx = 0:Lseg-1;
                for ip = 1:nP
                    P = P_list(ip);
                    lv = base - floor(idx / P);
                    c = dp(kprev, jc) + sum((ideal - double(lv)).^2);
                    if c < dp(kr, ic)
                        dp(kr, ic) = c;
                        prev(kr, ic) = j;
                        prevP(kr, ic) = P;
                        endL(kr, ic) = lv(end);
                    end
                end
            end
        end
    end

    costs = dp(2:K+1, N+1).';
    valid = costs < INF;
    if ~any(valid)
        error('nt2s DP failed to find a valid segmentation');
    end
    bestc = min(costs(valid));

    if force_nseg > 0
        nseg_used = min(max(1, force_nseg), K);
    else
        % smallest k capturing >=95% of SSE improvement toward best
        nseg_used = K;
        c1 = costs(1);
        if c1 > bestc + 1e-12
            for k = 1:K
                gain = (c1 - costs(k)) / (c1 - bestc);
                if gain >= 0.95
                    nseg_used = k;
                    break;
                end
            end
        else
            nseg_used = 1;
        end
    end

    segs = zeros(nseg_used, 2);
    k = nseg_used;
    i = N;
    while k > 0 && i > 0
        kr = k + 1;
        ic = i + 1;
        j = prev(kr, ic);
        P = prevP(kr, ic);
        segs(k, :) = [j + 1, P];
        i = j;
        k = k - 1;
    end
    seg_table = segs;
    best_lv = build_level_segs(N, L0, seg_table);
end
