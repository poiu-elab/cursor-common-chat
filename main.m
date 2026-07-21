clear all;
c       = physconst('LightSpeed');

N_chirp = 512;
fc0     = 75.9474e9;
tr0     = (48*54) * (1/54e6);  % 16us
deltaB  = 1.364501953e6;

v       = 10;
R       = 0;

nX0     = db2mag(20);
nAmp    = db2mag(-57);

tScale  = 1/54e6;
% tScale  = 1/756e6;

for v=0:10:100
    for n=1:N_chirp
        t3(n) = fc0*tr0/(fc0+(n-1)*deltaB);
        t2(n) = round(t3(n)/tScale)*tScale;
        t2(n) = tr0;
        st3(n)= n * fc0*tr0/(fc0+(n-1)*deltaB);
        if(n>1)
            nt3(n)= st3(n)-st3(n-1);
        elseif(n==1)
            nt3(n)= tr0;
        end
        nt2(n) = round(nt3(n)/tScale)*tScale;

        x0(n) = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(t2(1:n))) / c);
        x1(n) = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(t3(1:n))) / c);
        nx0(n) = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(nt2(1:n))) / c);
        nx1(n) = nX0 * exp(1j * 2 * pi * (fc0 + (n-1) * deltaB) * (2 * R + 2 * v * sum(nt3(1:n))) / c);
    end
    y0 = x0 + nAmp*randn(size(x0));
    y1 = x1 + nAmp*randn(size(x1));
    ny0 = nx0 + nAmp*randn(size(x0));
    ny1 = nx1 + nAmp*randn(size(x1));
    figure(1);
    subplot(121);
%     plot(t2);hold on;
%     plot(t3);
    plot(nt2);hold on;
    plot(nt3);
    hold off;
    title([sprintf('Tr0=%dus,N-chirp=%d',round(tr0*1e6),N_chirp)]);
    pause(0.01);
    subplot(122);
%     plot(mag2db(abs(fft(y0.*hanning(N_chirp).',N_chirp*1))));hold on;
%     plot(mag2db(abs(fft(y1.*hanning(N_chirp).',N_chirp*1))),'-o');hold on;
    plot(mag2db(abs(fft(ny0.*chebwin(N_chirp).',N_chirp*1))));hold on;
    plot(mag2db(abs(fft(ny1.*chebwin(N_chirp).',N_chirp*1))),'-o');hold on;
    hold off;
    grid on;grid minor;
    ylim([-40 80]);
    title([sprintf('scale=%dMHz,vel=%dm/s',1/tScale/1e6,v)]);
    pause(0.3);
end

%% Detect plateau lengths of nt2 (constant-value stair steps) and approximate period
% nt2 is quantized to tScale; a jump occurs when the quantized level changes.
dx = diff(nt2);
jump_indices = find(abs(dx) > tScale/2);   % last index of each plateau (except final)
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

% Approximate period from plateau lengths (exclude first/last which may be truncated)
inner_lengths = detected_lengths;
if numel(inner_lengths) >= 3
    inner_lengths = detected_lengths(2:end-1);
end
period_mean   = mean(inner_lengths);
period_median = median(inner_lengths);
% mode via histogram (works without Statistics Toolbox)
[counts, edges] = histcounts(inner_lengths, min(inner_lengths):max(inner_lengths)+1);
[~, imode] = max(counts);
period_mode = edges(imode);
% theory: local step length ≈ tScale / mean(|d nt3|)
period_theory = tScale / mean(abs(diff(nt3)));
approx_period = round(period_median);

fprintf('\n=== Approximate period (chirps per stair step) ===\n');
fprintf('mean(inner)   = %.3f\n', period_mean);
fprintf('median(inner) = %.3f\n', period_median);
fprintf('mode(inner)   = %d\n', period_mode);
fprintf('theory        = %.3f  (tScale/mean|d nt3|)\n', period_theory);
fprintf('approx_period = %d\n', approx_period);
fprintf('length sequence: %s\n', mat2str(detected_lengths));

figure(2);
subplot(211);
stem(detected_lengths, 'filled');
yline(approx_period, 'r--', sprintf('period≈%d', approx_period));
xlabel('plateau index'); ylabel('length (chirps)');
title(sprintf('nt2 constant-value lengths, approx period = %d', approx_period));
grid on;

subplot(212);
stairs(1:N_chirp, nt2*1e6);
xlabel('chirp index n'); ylabel('nt2 (us)');
title('nt2 staircase');
grid on;

