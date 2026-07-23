function [seg_table, best_lv, nseg_used, costs] = fit_nt2s_dp(levels, P_list, Nseg_max, force_nseg, fit_stride)
%FIT_NT2S_DP Piecewise integer-period fit minimizing SSE vs ideal levels.
%   DP chooses up to Nseg_max segments, each with one period from P_list.
%   force_nseg: 0 = auto elbow (>=95% SSE gain); >0 = exact segment count.
%   fit_stride: candidate boundary stride (1=exact; 2/4=faster approx).

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
            if isempty(j_list) || j_list(end) ~= i-1
                j_list = [j_list, i-1]; %#ok<AGROW>
            end
            for jj = 1:numel(j_list)
                j = j_list(jj);
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
