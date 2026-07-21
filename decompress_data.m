function dCmpData = decompress_data(DataIn,Params,DataOutputBitWidth,DataInputBitWidth,BlockHeadBitWidth,GroupHeadBitWidth) 

block_head = DataIn(1:BlockHeadBitWidth);
if(strcmp(block_head(1:2),'00')==1)
    Params.compressRate = 8/16;
elseif(strcmp(block_head(1:2),'01')==1)
    Params.compressRate = 7/16;
elseif(strcmp(block_head(1:2),'10')==1)
    Params.compressRate = 6/16;
elseif(strcmp(block_head(1:2),'11')==1)
    Params.compressRate = 5/16;
end
if(strcmp(block_head(3),'0')==1)
    Params.cp_mode_sel = 0;
else
    Params.cp_mode_sel = 1;
end
if(strcmp(block_head(4),'0')==1)
    Params.reshape_en = 0;
else
    Params.reshape_en = 1;
end

if(Params.reshape_en==0) 
    GroupDataBitWidth = floor(32*Params.antennas*Params.rangebins*Params.chirps*Params.compressRate-BlockHeadBitWidth-GroupHeadBitWidth); 
else
    GroupDataBitWidth = floor((32*Params.antennas*Params.rangebins*Params.chirps*Params.compressRate-BlockHeadBitWidth)/Params.rangebins-GroupHeadBitWidth);
end
dCmpData = zeros(Params.antennas*Params.rangebins,1);

if(Params.cp_mode_sel==0)
    dCmpData = Bfp_decompress(DataIn,Params,DataOutputBitWidth,DataInputBitWidth,BlockHeadBitWidth,GroupHeadBitWidth,GroupDataBitWidth);
else
	dCmpData = Ege_decompress(DataIn,Params,DataOutputBitWidth,DataInputBitWidth,BlockHeadBitWidth,GroupHeadBitWidth,GroupDataBitWidth);
end
dCmpData = complex(dCmpData(2:2:end),dCmpData(1:2:end));

if(Params.reshape_en==1)
    dCmpData = reshape(dCmpData,Params.antennas,Params.rangebins);
    dCmpData = dCmpData.';
    dCmpData = dCmpData(:);
end

end