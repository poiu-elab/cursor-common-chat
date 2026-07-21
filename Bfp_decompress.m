function DataOut = Bfp_decompress(Input,Params,DataOutputBitWidth,DataInputBitWidth,BlockHeadBitWidth,GroupHeadBitWidth,GroupDataBitWidth)

% Work with a logical bit vector and advance an index instead of repeatedly
% deleting Input(1:n), which copies the remaining character array each time.
if ischar(Input)
    bits = Input(:).' == '1';
elseif isa(Input, 'string')
    bits = char(Input);
    bits = bits(:).' == '1';
else
    bits = logical(Input(:).');
end

numSamples = Params.antennas * Params.rangebins;
DataOut = zeros(1, 2 * numSamples);
position = BlockHeadBitWidth + 1;

if Params.reshape_en == 0
    lsbI = nibble_value(bits, position);
    lsbQ = nibble_value(bits, position + 4);
    position = position + GroupHeadBitWidth;
    payloadWidth = floor(GroupDataBitWidth / (2 * numSamples));
    DataOut(:) = decode_components( ...
        bits, position, payloadWidth, numSamples, lsbI, lsbQ, ...
        DataOutputBitWidth / 2);
else
    payloadWidth = floor(GroupDataBitWidth / (2 * Params.antennas));
    valuesPerGroup = 2 * Params.antennas;
    for rangeIdx = 1:Params.rangebins
        lsbI = nibble_value(bits, position);
        lsbQ = nibble_value(bits, position + 4);
        position = position + GroupHeadBitWidth;
        outputIndices = (rangeIdx - 1) * valuesPerGroup ...
            + (1:valuesPerGroup);
        DataOut(outputIndices) = decode_components( ...
            bits, position, payloadWidth, Params.antennas, lsbI, lsbQ, ...
            DataOutputBitWidth / 2);
        position = position + valuesPerGroup * payloadWidth;
    end
end

end

function value = nibble_value(bits, firstBit)
value = 8 * double(bits(firstBit)) ...
      + 4 * double(bits(firstBit + 1)) ...
      + 2 * double(bits(firstBit + 2)) ...
      +     double(bits(firstBit + 3));
end

function values = decode_components( ...
    bits, firstBit, payloadWidth, numSamples, lsbI, lsbQ, componentWidth)

numComponents = 2 * numSamples;
lastBit = firstBit + numComponents * payloadWidth - 1;
payload = reshape(bits(firstBit:lastBit), payloadWidth, numComponents);
weights = 2 .^ (payloadWidth - 1:-1:0);
raw = weights * double(payload);

lsb = repmat([lsbI, lsbQ], 1, numSamples);
values = raw .* (2 .^ lsb);
negative = raw >= 2 ^ (payloadWidth - 1);
storedWidth = payloadWidth + lsb;
values(negative) = values(negative) ...
    + (2 .^ (componentWidth - storedWidth(negative)) - 1) ...
    .* (2 .^ storedWidth(negative));
wrap = values >= 2 ^ (componentWidth - 1);
values(wrap) = values(wrap) - 2 ^ componentWidth;
end


% function DataOut = Bfp_decompress(Input, Params, DataOutputBitWidth, DataInputBitWidth, BlockHeadBitWidth, GroupHeadBitWidth, GroupDataBitWidth)
% %BFP_DECOMPRESS_FAST 快速版本的BFP解压缩函数
% %   使用MEX文件实现，性能比原始MATLAB版本快10-50倍
%
% % 检查MEX文件是否存在，如果不存在则编译
% if ~exist('Bfp_decompress_mex', 'file')
%     compileMex();
% end
%
% % 确保输入数据格式正确
% if ischar(Input) || isstring(Input)
%     % 将字符串转换为uint8数组
%     Input = uint8(Input) - uint8('0');
% elseif islogical(Input)
%     % 逻辑数组直接使用
% else
%     % 尝试转换为uint8
%     Input = uint8(Input);
% end
%
% % 调用MEX函数
% DataOut = Bfp_decompress_mex(Input, Params, DataOutputBitWidth, DataInputBitWidth, BlockHeadBitWidth, GroupHeadBitWidth, GroupDataBitWidth);
%
% end
%
% function compileMex()
% %COMPILEMEX 编译MEX文件
%     fprintf('正在编译Bfp_decompress_mex...\n');
%
%     mexCmd = 'mex ';
%
%     % 根据平台设置编译选项
%     if ispc
%         mexCmd = [mexCmd, 'COMPFLAGS="$COMPFLAGS /O2" '];
%     else
%         mexCmd = [mexCmd, 'CXXOPTIMFLAGS="-O2" '];
%     end
%
%     mexCmd = [mexCmd, 'Bfp_decompress_mex.cpp'];
%
%     try
%         eval(mexCmd);
%         fprintf('编译成功！\n');
%     catch ME
%         error('编译失败: %s\n请确保已安装MATLAB编译器并配置了C++编译器。', ME.message);
%     end
% end
