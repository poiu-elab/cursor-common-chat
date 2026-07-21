function rfft_decomp_result = rfft_compress(data,nSample,nChirp,nRx,nTx,Params)
    blockSize = Params.rangebins * Params.chirps * Params.antennas;
    blockBytes = blockSize * 2;
    stepLength = (128 - (mod(blockBytes * 8 - 1, 128) + 1)) / 8;
    blockStride = blockBytes + stepLength;

    if mod(nSample, Params.rangebins) ~= 0 ...
            || mod(nChirp, Params.chirps) ~= 0 ...
            || mod(nRx, Params.antennas) ~= 0
        error('rfft_compress:InvalidDimensions', ...
            'Input dimensions must be divisible by the BFP block dimensions.');
    end

    numRangeBlocks = nSample / Params.rangebins;
    numChirpBlocks = nChirp / Params.chirps;
    numAntennaBlocks = nRx / Params.antennas;
    numBlocks = nTx * numRangeBlocks * numChirpBlocks * numAntennaBlocks;
    requiredBytes = (numBlocks - 1) * blockStride + blockBytes;
    if numel(data) < requiredBytes
        error('rfft_compress:InsufficientData', ...
            'Expected at least %d input bytes, but received %d.', ...
            requiredBytes, numel(data));
    end

    data = uint8(data(:));
    cubeSize = [nSample, nChirp, nRx, nTx];
    DcmpData = complex(zeros(cubeSize));
    encode_data = complex(zeros(cubeSize));
    [localRange, localChirp, localAntenna] = ind2sub( ...
        [Params.rangebins, Params.chirps, Params.antennas], ...
        (1:blockSize).');

    % Limit temporary bit matrices while still decoding thousands of blocks
    % per call. The old path called dec2bin/bin2dec once for every block.
    targetTemporaryBits = 2 ^ 21;
    blocksPerChunk = max(1, floor(targetTemporaryBits / (blockBytes * 8)));
    for firstBlock = 1:blocksPerChunk:numBlocks
        lastBlock = min(firstBlock + blocksPerChunk - 1, numBlocks);
        blockNumbers = firstBlock:lastBlock;
        blockStarts = 1 + (blockNumbers - 1) * blockStride;
        byteIndices = bsxfun(@plus, (0:blockBytes-1).', blockStarts);
        rawBlocks = data(byteIndices);
        decodedChunk = complex(zeros(blockSize, numel(blockNumbers)));

        % cp_mode_sel is the third MSB of the first byte after word reversal.
        isEge = bitget(rawBlocks(8, :), 6) ~= 0;
        bfpColumns = find(~isEge);
        if ~isempty(bfpColumns)
            decodedChunk(:, bfpColumns) = bfp_decompress_blocks( ...
                rawBlocks(:, bfpColumns), ...
                Params, Params.DataOutputBitWidth_Dcm);
        end

        % Preserve the legacy EGE dispatch when that external decoder exists.
        % BFP data, which is the normal path, never enters this scalar loop.
        egeColumns = find(isEge);
        for localColumn = reshape(egeColumns, 1, [])
            orderedBytes = reshape(rawBlocks(:, localColumn), 8, []);
            orderedBytes = flip(orderedBytes, 1);
            bitCharacters = dec2bin(double(orderedBytes(:)), 8);
            encodedBits = reshape(bitCharacters.', 1, []);
            decodedChunk(:, localColumn) = decompress_data( ...
                encodedBits, Params, Params.DataOutputBitWidth_Dcm, ...
                Params.DataInputBitWidth_Dcm, Params.BlockHeadBitWidth, ...
                Params.GroupHeadBitWidth);
        end

        encodedChunk = complex( ...
            double(rawBlocks(1:2:end, :)), ...
            double(rawBlocks(2:2:end, :)));

        [antennaBlock, chirpBlock, rangeBlock, txBlock] = ind2sub( ...
            [numAntennaBlocks, numChirpBlocks, numRangeBlocks, nTx], ...
            blockNumbers);
        globalRange = bsxfun(@plus, localRange, ...
            (rangeBlock - 1) * Params.rangebins);
        globalChirp = bsxfun(@plus, localChirp, ...
            (chirpBlock - 1) * Params.chirps);
        globalAntenna = bsxfun(@plus, localAntenna, ...
            (antennaBlock - 1) * Params.antennas);
        linearIndices = globalRange ...
            + (globalChirp - 1) * nSample ...
            + (globalAntenna - 1) * nSample * nChirp;
        linearIndices = bsxfun(@plus, linearIndices, ...
            (txBlock - 1) * nSample * nChirp * nRx);
        DcmpData(linearIndices) = decodedChunk;
        encode_data(linearIndices) = encodedChunk;
    end

    % Keep the original API: when nTx > 1, the legacy loop returned the
    % final transmitter because each transmitter overwrote the same cube.
    DcmpData = DcmpData(:, :, :, end);
    encode_data = encode_data(:, :, :, end);

    rfft_decomp_result.input = encode_data;
    rfft_decomp_result.output = DcmpData;
end