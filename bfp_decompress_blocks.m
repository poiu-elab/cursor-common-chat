function DataOut = bfp_decompress_blocks(rawBlocks, Params, DataOutputBitWidth)
%BFP_DECOMPRESS_BLOCKS Decode multiple byte-packed BFP blocks at once.
%   RAWBLOCKS is a byte matrix with one compressed block per column. Bytes
%   are reordered in groups of eight to match the capture stream format.

if Params.chirps ~= 1
    error('bfp_decompress_blocks:UnsupportedChirpBlock', ...
        'The existing BFP format supports Params.chirps == 1 only.');
end
if Params.BlockHeadBitWidth < 4 || Params.GroupHeadBitWidth < 8
    error('bfp_decompress_blocks:InvalidHeaderWidth', ...
        'BFP block and group headers require at least 4 and 8 bits.');
end
if mod(size(rawBlocks, 1), 8) ~= 0
    error('bfp_decompress_blocks:InvalidBlockSize', ...
        'The compressed block size must be a multiple of eight bytes.');
end
if mod(DataOutputBitWidth, 2) ~= 0
    error('bfp_decompress_blocks:InvalidOutputWidth', ...
        'DataOutputBitWidth must be even.');
end

numBlocks = size(rawBlocks, 2);
numSamples = Params.antennas * Params.rangebins;
bytesPerBlock = size(rawBlocks, 1);

% The capture stream stores each eight-byte word in reverse byte order.
orderedBytes = reshape(uint8(rawBlocks), 8, [], numBlocks);
orderedBytes = flip(orderedBytes, 1);
orderedBytes = reshape(orderedBytes, bytesPerBlock, numBlocks);

% Expand all bytes in a chunk together. This replaces dec2bin and repeated
% character concatenation in the original per-block implementation.
bits = false(bytesPerBlock * 8, numBlocks);
for bitIdx = 1:8
    bits(bitIdx:8:end, :) = bitget(orderedBytes, 9 - bitIdx) ~= 0;
end

compressionCode = 2 * double(bits(1, :)) + double(bits(2, :));
if any(bits(3, :))
    error('bfp_decompress_blocks:UnsupportedMode', ...
        'Only BFP blocks (cp_mode_sel == 0) can be batch decoded.');
end
reshapeMode = double(bits(4, :));

DataOut = complex(zeros(numSamples, numBlocks));
compressionRates = [8, 7, 6, 5] / 16;

for mode = 0:1
    for code = 0:3
        columns = find(reshapeMode == mode & compressionCode == code);
        if isempty(columns)
            continue;
        end

        rate = compressionRates(code + 1);
        if mode == 0
            groupDataBitWidth = floor(32 * numSamples * rate ...
                - Params.BlockHeadBitWidth - Params.GroupHeadBitWidth);
            payloadWidth = floor(groupDataBitWidth / (2 * numSamples));

            groupStart = Params.BlockHeadBitWidth + 1;
            lsbI = nibble_value(bits, groupStart, columns);
            lsbQ = nibble_value(bits, groupStart + 4, columns);
            DataOut(:, columns) = decode_payload( ...
                bits, groupStart + Params.GroupHeadBitWidth, ...
                payloadWidth, numSamples, columns, lsbI, lsbQ, ...
                DataOutputBitWidth / 2);
        else
            groupDataBitWidth = floor( ...
                (32 * numSamples * rate - Params.BlockHeadBitWidth) ...
                / Params.rangebins - Params.GroupHeadBitWidth);
            payloadWidth = floor(groupDataBitWidth / (2 * Params.antennas));
            groupStride = Params.GroupHeadBitWidth ...
                + 2 * Params.antennas * payloadWidth;

            grouped = complex(zeros( ...
                Params.antennas, Params.rangebins, numel(columns)));
            for rangeIdx = 1:Params.rangebins
                headerStart = Params.BlockHeadBitWidth + 1 ...
                    + (rangeIdx - 1) * groupStride;
                lsbI = nibble_value(bits, headerStart, columns);
                lsbQ = nibble_value(bits, headerStart + 4, columns);
                values = decode_payload(bits, ...
                    headerStart + Params.GroupHeadBitWidth, ...
                    payloadWidth, Params.antennas, columns, lsbI, lsbQ, ...
                    DataOutputBitWidth / 2);
                grouped(:, rangeIdx, :) = reshape(values, ...
                    Params.antennas, 1, numel(columns));
            end

            % Match decompress_data's antenna/range transpose for reshape_en.
            DataOut(:, columns) = reshape( ...
                permute(grouped, [2, 1, 3]), numSamples, numel(columns));
        end
    end
end

end

function value = nibble_value(bits, firstBit, columns)
value = 8 * double(bits(firstBit, columns)) ...
      + 4 * double(bits(firstBit + 1, columns)) ...
      + 2 * double(bits(firstBit + 2, columns)) ...
      +     double(bits(firstBit + 3, columns));
end

function values = decode_payload( ...
    bits, firstBit, payloadWidth, numSamples, columns, lsbI, lsbQ, ...
    componentWidth)

numColumns = numel(columns);
lastBit = firstBit + 2 * numSamples * payloadWidth - 1;
payload = bits(firstBit:lastBit, columns);
payload = reshape(payload, payloadWidth, 2 * numSamples, numColumns);
weights = reshape(2 .^ (payloadWidth - 1:-1:0), [], 1, 1);
raw = reshape(sum(bsxfun(@times, double(payload), weights), 1), ...
    2 * numSamples, numColumns);

lsb = zeros(2 * numSamples, numColumns);
lsb(1:2:end, :) = repmat(lsbI, numSamples, 1);
lsb(2:2:end, :) = repmat(lsbQ, numSamples, 1);

% Preserve the original fixed-width sign extension, including wraparound
% behavior for malformed headers whose mantissa plus LSB exceeds 16 bits.
signedValue = raw .* (2 .^ lsb);
negative = raw >= 2 ^ (payloadWidth - 1);
storedWidth = payloadWidth + lsb;
signedValue(negative) = signedValue(negative) ...
    + (2 .^ (componentWidth - storedWidth(negative)) - 1) ...
    .* (2 .^ storedWidth(negative));
wrap = signedValue >= 2 ^ (componentWidth - 1);
signedValue(wrap) = signedValue(wrap) - 2 ^ componentWidth;

% Stored component order is imaginary, real.
values = complex(signedValue(2:2:end, :), signedValue(1:2:end, :));
end
