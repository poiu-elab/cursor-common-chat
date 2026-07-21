function DataOut = Bfp_decompress(Input,Params,DataOutputBitWidth,DataInputBitWidth,BlockHeadBitWidth,GroupHeadBitWidth,GroupDataBitWidth)

Input(1:BlockHeadBitWidth) = [];
group_head                 = Input(1:GroupHeadBitWidth);
lsbI                       = bin2dec(group_head(1:4));
lsbQ                       = bin2dec(group_head(5:8));
Input(1:GroupHeadBitWidth) = [];
DataOut                    = zeros(1,Params.antennas*Params.rangebins*2);
cnt                        = 1;

if(Params.reshape_en==0)  
    lenI = lsbI+floor(GroupDataBitWidth/(2*Params.antennas*Params.rangebins));
    lenQ = lsbQ+floor(GroupDataBitWidth/(2*Params.antennas*Params.rangebins));
    for ind = 1:Params.antennas*Params.rangebins
        % imag part
        data               = bin2dec(Input(1:lenI-lsbI));
        Input(1:lenI-lsbI) = [];
        if(data<2^(lenI-lsbI-1))
            data  = data*2^lsbI;
        else
            data  = data*2^lsbI+(2^(DataOutputBitWidth/2-lenI)-1)*2^lenI;
        end
        if(data>=2^(DataOutputBitWidth/2-1))
            data = data-2^(DataOutputBitWidth/2);
        end
        DataOut(cnt) = data;
        cnt          = cnt+1;
        % real part
        data               = bin2dec(Input(1:lenQ-lsbQ));
        Input(1:lenQ-lsbQ) = [];
        if(data<2^(lenQ-lsbQ-1))
            data  = data*2^lsbQ;
        else
            data  = data*2^lsbQ+(2^(DataOutputBitWidth/2-lenQ)-1)*2^lenQ;
        end
        if(data>=2^(DataOutputBitWidth/2-1))
            data = data-2^(DataOutputBitWidth/2);
        end
        DataOut(cnt) = data;
        cnt          = cnt+1;
    end
else
    for BcntIdx = 1:Params.rangebins
        lenI = lsbI + floor(GroupDataBitWidth/2/Params.antennas);
        lenQ = lsbQ + floor(GroupDataBitWidth/2/Params.antennas);
        for ind = 1:Params.antennas 
            % imag part
            data               = bin2dec(Input(1:lenI-lsbI));
            Input(1:lenI-lsbI) = [];
            if(data<2^(lenI-lsbI-1))
                data  = data*2^lsbI;
            else
                data  = data*2^lsbI+(2^(DataOutputBitWidth/2-lenI)-1)*2^lenI;
            end
            if(data>=2^(DataOutputBitWidth/2-1))
                data = data-2^(DataOutputBitWidth/2);
            end
            DataOut(cnt) = data;
            cnt          = cnt+1;
            % real part
            data               = bin2dec(Input(1:lenQ-lsbQ));
            Input(1:lenQ-lsbQ) = [];
            if(data<2^(lenQ-lsbQ-1))
                data  = data*2^lsbQ;
            else
                data  = data*2^lsbQ+(2^(DataOutputBitWidth/2-lenQ)-1)*2^lenQ;
            end
            if(data>=2^(DataOutputBitWidth/2-1))
                data = data-2^(DataOutputBitWidth/2);
            end
            DataOut(cnt) = data;
            cnt          = cnt+1;
        end
        if BcntIdx < Params.rangebins
            group_head                 = Input(1:GroupHeadBitWidth);
            lsbI                       = bin2dec(group_head(1:4));
            lsbQ                       = bin2dec(group_head(5:8));
            Input(1:GroupHeadBitWidth) = [];
        end
    end 
end

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
