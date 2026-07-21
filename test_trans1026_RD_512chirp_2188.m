clear all;
isSingleFile  = 1;  %% 1：单个文�?                       0：文件夹
saveCube      = 0;  %% 1：存储解压后数据                  0：不存储解压后数�?
RDdisplay     = 1;  %% 1：显示RD�?                       0：不显示
isSingleFrame = 1;  %% 1：单帧数�?                       0：多帧数�?
saveCoef      = 0;  %% 1：存储系数文件（生成幅相系数�?   0：不使能（直接读取已有的幅相系数�?
angleAnalyse  = 1;  %% 1：分析角度�?�?                   0：不分析
bpmOn         = 0;  %% 1：bpm使能                         0：不使能
boardType     = 30125;     %%30125：斯凯瑞丽采集板


nCplxSample = 496;
nRx = 8;

%[FileName,PathName] = uigetfile('8088-512-LKVCO-oldsdk-880486f5.txt','select data file'); 
%
%
%if isSingleFile == 1
%    FilePath = strcat(PathName,FileName);
%    filecnt = 1;
%else
%    Filename = dir([PathName '*.dat']);
%    filecnt = length(Filename);
%end
tic;
%FilePath = '.\nf\8088-768-HKVCO-oldsdk-7c7099e93_new.bin';  nChirp = 768;
%FilePath = '.\nf\8088-768-LKVCO-oldsdk-edbfcb08.bin';       nChirp = 768;
%FilePath = '.\nf\8088-512-LKVCO-oldsdk-880486f5.bin';       nChirp = 512;
%FilePath = '.\100m\8088-768-HKVCO-oldsdk-7c7099e93.bin';    nChirp = 768;
%FilePath = '.\100m\8088-768-LKVCO-oldsdk-edbfcb08.bin';     nChirp = 768;
%FilePath = '.\100m\8088-512-LKVCO-oldsdk-880486f5.bin';     nChirp = 512;
FilePath = '.\200m\8088-768-HKVCO-oldsdk-7c7099e93.bin';    nChirp = 768;
%FilePath = '.\200m\8088-768-LKVCO-oldsdk-edbfcb08.bin';     nChirp = 768;
%FilePath = '.\200m\8088-512-LKVCO-oldsdk-880486f5.bin';     nChirp = 512;

for fileNum=1:1
    if isSingleFile ~= 1
        FilePath = strcat(PathName,Filename(fileNum).name);
    end
    
    fid = fopen(FilePath,'rb');
    dataTest = fread(fid, 'uint32');
    fclose(fid);
    
    totalLen = nCplxSample*nChirp*nRx*4/2 + 32 + 32;    %datalength + frame_head(32byte) + frame_end(32byte)
    dataPos = find(dataTest == hex2dec('AAAA5555'));
    if isSingleFrame == 1
        dataPos(2) = 100;
        fid = fopen(FilePath,'rb');
        fseek(fid,0,'bof');
        frameNum = 1;
    end

    for nn  = 1:frameNum
        dataHead1 = dec2hex(fread(fid,16,'uint8'));
        frameCntALL(:,nn) = fread(fid,1,'uint32');
        CarSpeedALL(:,nn) = fread(fid,1,'float');
        reserve = fread(fid,8,'uint8');
        dataCompressNarrowALL(:,nn) = fread(fid,nChirp*nCplxSample*8*2, 'uint8');
        dataEnd =  dec2hex(fread(fid,32,'uint8'));
    end
    fclose(fid);

    for cnt = 1:frameNum
        Params.rangebins              = 4; 
        Params.chirps                 = 1; 
        Params.antennas               = 8;
        Params.DataOutputBitWidth_Dcm = 32; 
        Params.DataInputBitWidth_Dcm  = 128;  
        Params.BlockHeadBitWidth      = 4; 
        Params.GroupHeadBitWidth      = 8;
        data = dataCompressNarrowALL(:,cnt);
        rfft_decomp_result            = rfft_compress(data,nCplxSample,nChirp,nRx,1,Params);
        rfftCube = rfft_decomp_result.output;
        rfftCubeTrans = rfftCube;
        if boardType==30125
            rfftCube(:,:,3) = rfftCubeTrans(:,:,5);
            rfftCube(:,:,4) = rfftCubeTrans(:,:,6);
            rfftCube(:,:,5) = rfftCubeTrans(:,:,3);
            rfftCube(:,:,6) = rfftCubeTrans(:,:,4);
        end
        save([sprintf('rfftCube-%s.mat',char(datetime('now', 'Format', 'yyyyMMdd-HHmmss')))],'rfftCube');
        commonCfg.useRange  =   size(rfftCube,1);
        commonCfg.numChirp  =   size(rfftCube,2);
        commonCfg.numRx     =   size(rfftCube,3);
        abuf_fft2d          = fft(rfftCube .* repmat(hanning(commonCfg.numChirp).', commonCfg.useRange, 1, commonCfg.numRx),commonCfg.numChirp*4,2);
        abuf_fft2d_abs_mean = squeeze(mean(abs(abuf_fft2d),2));
        
        figure(100);
        meshz(mag2db(abs(abuf_fft2d(:,:,1))));
        view([90 0]);
        figure(101);
        plot(400/commonCfg.useRange*[0:commonCfg.useRange-1],mag2db(abuf_fft2d_abs_mean)); %ylim([50 90]); xlim([0 400]);
        hold on;
        hold off;ylim([25 60]);
        xlabel('����(m)');
        ylabel('noise floor');
        grid on;grid minor;

        d1Idx=133; d2Idx=1866;   % 27M 512
%        d1Idx=128; d2Idx=1539;   % 54M 768
        s=abs(abuf_fft2d(d1Idx,d2Idx,1)); n=mean(abs(abuf_fft2d(d1Idx,[1:d2Idx-20 d2Idx+20:end],1)));
        figure(103);
        plot(mag2db(abs(abuf_fft2d(d1Idx,:,1))));
        title([sprintf('r=%gm \n snr=%g',d1Idx*350/commonCfg.useRange,mag2db(s/n))]);
        grid on;grid minor;

        d1Idx=265; d2Idx=1683;   % 27M 512
%        d1Idx=254; d2Idx=16;   % 54M 768
        s=abs(abuf_fft2d(d1Idx,d2Idx,1)); n=mean(abs(abuf_fft2d(d1Idx,[1:d2Idx-20 d2Idx+20:end],1)));
        figure(104);
        plot(mag2db(abs(abuf_fft2d(d1Idx,:,1))));
        title([sprintf('r=%gm \n snr=%g',d1Idx*350/commonCfg.useRange,mag2db(s/n))]);
        grid on;grid minor;

%        if saveCube == 1
%            saveFilePath = strcat(PathName,strrep(FilePath,'.dat',''),num2str(cnt),'_',num2str(frameCntALL(cnt)),'_',num2str(CarSpeedALL(cnt)),'.mat');
%            save(saveFilePath,'rfftCube');
%        end
%        
%        if RDdisplay==1
%            cubeSize               = [nCplxSample, nChirp, 8];
%            d2Cube = zeros(cubeSize(1),cubeSize(2),cubeSize(3));
%            d2Cubetdm = zeros(cubeSize(1),cubeSize(2)/8,8,cubeSize(3));
%            
%            if bpmOn == 1
%                load('bpmcode.mat');
%                for cnttt = 1:60
%                    tmp = squeeze(rfftCube(1,:,1));
%                    d2Win = bpmcode(cnttt,:).*hanning(cubeSize(2)).';
%                    dataFFT2DSel(cnttt) = max(max(abs(fft(tmp.*d2Win))));
%    %                 dataFFT2DSel(cnttt) = max(max(abs(fft2(dataChk(1:nCplxSample,:).*repmat(hanning(nCplxSample),1,nChirp).*repmat((bpmcode(cnttt,:).*hanning(nChirp).'),nCplxSample,1),512,768))));
%                end
%                [val,indexbpm] = max(dataFFT2DSel);             
%                d2Win = bpmcode(indexbpm,:).*hanning(cubeSize(2)).';
%            else
%                %d2Win = hanning(cubeSize(2)).';
%                d2Win = hanning(cubeSize(2)/8).';%%tdm 64 chirp
%            end
%            for i = 1: cubeSize(1)         %range
%                for j = 1: 8               % TX
%                    for k = 1:cubeSize(3)  %RX
%                        tmp = squeeze(rfftCube(i,1+64*(j-1):64*j,k));
%                        %tmp =squeeze(rfftCube(i,j:8:end,k));
%                        fftOut = fft(tmp.*d2Win);
%                        %d2Cube(i,:,j) = fftOut;
%                        d2Cubetdm(i,:,j,k) = fftOut;
%                    end
%                end
%            end
%            figure(1);   
%            for j = 1:8
%                subplot(2,4,j);
%                modsum = sum(abs(d2Cubetdm(:,:,j,:)),4);
%                %imagesc(db(modsum).');
%                mesh(db(modsum));
%                title(['TX',num2str(j)]);
%            end
%            figure(2);
%            for k = 1:8
%                subplot(2,4,k);
%                j = 4;                %TX
%                mod = abs(d2Cubetdm(:,:,j,k));
%                %imagesc(db(modsum).');
%                mesh(db(mod));
%                title(['TX',num2str(j),'RX',num2str(k)]);
%            end
%            figure(3);
%            for k = 1:8
%                subplot(2,4,k);
%                j = 5;                %TX
%                mod = abs(d2Cubetdm(:,:,j,k));
%                %imagesc(db(modsum).');
%                mesh(db(mod));
%                title(['TX',num2str(j),'RX',num2str(k)]);
%            end
%        end
    end
%    Hist_noise = getHistNoise(modsum);
%    figure(222)
%    plot(db(Hist_noise))
%    grid on
%    grid minor
%    hold on
%    %plot(db(modsum(:,169)));
%    hold off
%    figure(3),mesh(db(modsum));
%    %% calc SNR
%    modsum_db = db(modsum);
%    [r_max,d_max]=find(modsum_db == max(max(modsum_db(100:end,:))));
%    SNR = modsum_db(r_max,d_max) - db(Hist_noise(r_max));
%    fprintf('r:%d d:%d fudu:%f SNR:%f noise:%f \n',r_max,d_max,modsum_db(r_max,d_max),SNR,db(Hist_noise(r_max)));
%    %% plot
%    figure(222),
%    hold on,
%    plot(db(modsum(:,d_max))),
%    title(FileName);
%    %%plot
%    figure(333),
%    plot(db(modsum(r_max,:)));
%    figure(100);
%    plot(unwrap(angle(rfftCube(133,:,2))),'--o');
end
toc;
