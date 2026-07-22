clear all;
c       = physconst('LightSpeed');

N_chirp = 512;
fc0     = 75.9474e9;
tr0     = (48*54) * (1/54e6);  % 48us
deltaB  = 1.364501953e6;

v       = 10;
R       = 0;

nX0     = db2mag(20);
nAmp    = db2mag(-57);

tScale  = 1/54e6;
% tScale  = 1/756e6;

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

%% Detect plateau lengths of nt2 and approximate a base period
dx = diff(nt2);
jump_indices = find(abs(dx) > tScale/2);
plateau_starts = [1, jump_indices + 1];
plateau_ends   = [jump_indices, N_chirp];
detected_lengths = plateau_ends - plateau_starts + 1;
plateau_values   = nt2(plateau_starts);

fprintf('\n=== nt2 plateau lengths (constant-value runs) ===\n');
fprintf('tScale = %.6f ns, plateaus = %d\n', tScale*1e9, numel(detected_lengths));
fprintf('%4s %6s %6s %6s %14s\n', 'idx', 'start', 'end', 'len', 'nt2(us)');
for k = 1:numel(detected_lengths)
    fprintf('%4d %6d %6d %6d %14.9f\n', k, plateau_starts(k), plateau_ends(k), ...
        detected_lengths(k), plateau_values(k)*1e6);
end

inner_lengths = detected_lengths;
if numel(inner_lengths) >= 3
    inner_lengths = detected_lengths(2:end-1);
end
period_mean   = mean(inner_lengths);
period_median = median(inner_lengths);
[counts, edges] = histcounts(inner_lengths, min(inner_lengths):max(inner_lengths)+1);
[~, imode] = max(counts);
period_mode = edges(imode);
period_theory = tScale / mean(abs(diff(nt3)));
approx_period = round(period_median);

fprintf('\n=== Approximate period (chirps per stair step) ===\n');
fprintf('mean(inner)   = %.3f\n', period_mean);
fprintf('median(inner) = %.3f\n', period_median);
fprintf('mode(inner)   = %d\n', period_mode);
fprintf('theory        = %.3f  (tScale/mean|d nt3|)\n', period_theory);
fprintf('approx_period = %d\n', approx_period);
fprintf('length sequence: %s\n', mat2str(detected_lengths));

%% Build nt2s: piecewise fixed-integer-period staircase (easy to implement)
% Each segment uses one whole period P: level drops by 1 every P chirps.
% Example (auto result for current params):
%   chirp 1..111  -> P=10
%   chirp 112..512 -> P=11
% Implement as: if n < seg_table(2,1), use P1, else use P2.
levels = round(nt2 / tScale);
L0 = levels(1);
Nseg_max = 2;                          % set 1 for single period only
P_list = unique([approx_period-1, approx_period, approx_period+1]);
P_list = P_list(P_list >= 1);

% Manual override (skip search): use_manual_segs=1
use_manual_segs = 0;
manual_segs = [1, 10; 112, 11];        % [start_chirp, period]

if use_manual_segs
    best_segs = manual_segs;
    best_lv = build_level_segs(N_chirp, L0, best_segs);
else
    best_err = inf;
    best_segs = [1, approx_period];
    best_lv = build_level_segs(N_chirp, L0, best_segs);
    nP = numel(P_list);
    for K = 1:Nseg_max
        if K == 1
            brk_list = 0;
        else
            brk_list = 2:(N_chirp-1);
        end
        for brk = brk_list
            nComb = nP^K;
            for ic = 0:nComb-1
                tmp = ic;
                Pseg = zeros(1, K);
                for ks = 1:K
                    Pseg(ks) = P_list(mod(tmp, nP) + 1);
                    tmp = floor(tmp / nP);
                end
                if K == 1
                    segs = [1, Pseg(1)];
                else
                    segs = [1, Pseg(1); brk, Pseg(2)];
                end
                lv = build_level_segs(N_chirp, L0, segs);
                err = sqrt(mean((double(levels) - double(lv)).^2));
                err = err + 0.05 * abs(lv(end) - levels(end));
                if err < best_err
                    best_err = err;
                    best_segs = segs;
                    best_lv = lv;
                end
            end
            if K == 1
                break;
            end
        end
    end
end

nt2s = best_lv * tScale;
seg_table = best_segs;

fprintf('\n=== nt2s segment table (start chirp, period) ===\n');
disp(seg_table);
fprintf('nt2s level RMSE vs nt2 = %.4f, end level %d (ideal %d)\n', ...
    sqrt(mean((double(levels)-double(best_lv)).^2)), best_lv(end), levels(end));
fprintf('cumsum final err = %.6f us\n', (sum(nt2s)-sum(nt2))*1e6);

% lengths of the regularized staircase
dxs = diff(nt2s);
jump_s = find(abs(dxs) > tScale/2);
len_s = diff([0, jump_s, N_chirp]);
fprintf('nt2s plateau lengths: %s\n', mat2str(len_s));

figure(2);
subplot(311);
stem(detected_lengths, 'filled'); hold on;
plot([1, numel(detected_lengths)], [approx_period, approx_period], 'r--');
hold off;
xlabel('plateau index'); ylabel('length');
title(sprintf('nt2 lengths (irregular), base period≈%d', approx_period));
grid on;

subplot(312);
stem(len_s, 'filled');
xlabel('plateau index'); ylabel('length');
title('nt2s lengths (segment-regularized whole periods)');
grid on;

subplot(313);
stairs(1:N_chirp, nt2*1e6, 'b'); hold on;
stairs(1:N_chirp, nt2s*1e6, 'r');
% mark segment boundaries
yl = [min(nt2s) min(nt2) max(nt2s) max(nt2)] * 1e6;
for iseg = 2:size(seg_table, 1)
    plot([seg_table(iseg, 1), seg_table(iseg, 1)], [min(yl), max(yl)], 'k--');
end
hold off;
xlabel('chirp index n'); ylabel('us');
legend('nt2 (ideal quant)', 'nt2s (regularized)', 'Location', 'best');
title('nt2 vs nt2s');
grid on;

%% Velocity sweep evaluation: nt2 vs nt2s vs continuous nt3
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
    plot(nt2); hold on;
    plot(nt2s); plot(nt3);
    hold off;
    legend('nt2', 'nt2s', 'nt3', 'Location', 'best');
    title(sprintf('Tr0=%dus,N-chirp=%d', round(tr0*1e6), N_chirp));

    subplot(122);
    w = chebwin(N_chirp).';
    plot(mag2db(abs(fft(ny0 .*w, N_chirp*1)))); hold on;
    plot(mag2db(abs(fft(ny0s.*w, N_chirp*1))), '-.');
    plot(mag2db(abs(fft(ny1 .*w, N_chirp*1))), '-o');
    hold off;
    legend('nt2', 'nt2s', 'nt3', 'Location', 'best');
    grid on; grid minor;
    ylim([-40 80]);
    title(sprintf('scale=%dMHz,vel=%dm/s', 1/tScale/1e6, v));
    pause(0.3);
end

%% Local helper: build level staircase from segment table
% segs: Mx2, each row [start_chirp_1based, period]
function lv = build_level_segs(N, L0, segs)
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
