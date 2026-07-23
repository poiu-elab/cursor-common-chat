function lv = build_level_segs(N, L0, segs)
%BUILD_LEVEL_SEGS Build staircase levels from segment table.
%   segs: Mx2, each row [start_chirp_1based, period]
%   Each segment resets the chirp counter; level drops by 1 every P chirps.

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
