clc; clear; close all;

fid = fopen("D:\ÁÙÊ±Êý¾Ý\1768481790323251088.bin", 'rb');
if (fid < 0)
    fprintf("Open bin file failed.\n");
end

% data = fread(fid, cnt, 'int16');
data = fread(fid,'int16');
fclose(fid);

data_q = data(1:2:end);
data_i = data(2:2:end);
rfftData_tmp_i = reshape(data_i,4,8,[]);
rfftData_tmp_q = reshape(data_q,4,8,[]);
rfftData = zeros(496,768,8);
for i = 1: 496/4
    for j = 1:768
        for k = 1:8
            rfftData(4*(i-1)+1:4*i,j,k) = complex(rfftData_tmp_i(:,k,(i-1)*768+j),rfftData_tmp_q(:,k,(i-1)*768+j));
        end
    end
end



cubeSize = size(rfftData);
d2Cube = zeros(cubeSize(1),cubeSize(2),cubeSize(3));
d2Win = hanning(cubeSize(2)).';
for i = 1: cubeSize(1)
    for j = 1: cubeSize(3)
        tmp = squeeze(rfftData(i,:,j));
        fftOut = fft(tmp.*d2Win);
        d2Cube(i,:,j) = fftOut;
    end
end
figure;
for i = 1:cubeSize(3)
   mesh(abs(squeeze(d2Cube(:,:,i))));
    hold on;
end
figure;
% for i = 1:cubeSize(3)
%    mesh(abs(squeeze(d2Cube(:,:,i))));
%     hold on;
% end
modsum = sum(d2Cube,3);
mesh(db(modsum))
