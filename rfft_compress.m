function rfft_decomp_result = rfft_compress(data,nSample,nChirp,nRx,nTx,Params)
    blockSize   = Params.rangebins*Params.chirps*Params.antennas;
    step_length = (128 - (mod(blockSize * 2 * 8 - 1, 128) + 1)) / 8;
    cnt         = 1;
    for n = 1:nTx
        for i = 1:Params.rangebins:nSample
            for j = 1:Params.chirps:nChirp
                for m = 1:Params.antennas:nRx
                    encodeData = [];
                    temp_data  = data(cnt:cnt+blockSize*2-1);
                    cnt        = cnt+step_length+blockSize*2;
                    for k = 1:8:blockSize*2
                        encodeData = [encodeData,dec2bin(temp_data(k+7),8),dec2bin(temp_data(k+6),8),dec2bin(temp_data(k+5),8),dec2bin(temp_data(k+4),8),...
                                                 dec2bin(temp_data(k+3),8),dec2bin(temp_data(k+2),8),dec2bin(temp_data(k+1),8),dec2bin(temp_data(k+0),8)];
                    end
                    dCmpData = decompress_data(encodeData,Params,Params.DataOutputBitWidth_Dcm,Params.DataInputBitWidth_Dcm,Params.BlockHeadBitWidth,Params.GroupHeadBitWidth); 
                    DcmpData(i:i+Params.rangebins-1,j:j+Params.chirps-1,m:m+Params.antennas-1)    = reshape(dCmpData,[Params.rangebins,Params.chirps,Params.antennas]);
                    temp_data_cplx                                                                = complex(temp_data(1:2:end),temp_data(2:2:end));
                    encode_data(i:i+Params.rangebins-1,j:j+Params.chirps-1,m:m+Params.antennas-1) = reshape(temp_data_cplx,[Params.rangebins,Params.chirps,Params.antennas]);
                end
            end
        end
    end
    rfft_decomp_result.input = encode_data;
    rfft_decomp_result.output = DcmpData;   
end
    
   