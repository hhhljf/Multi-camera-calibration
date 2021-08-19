function uvc = readvc(path)

% 读取单个.vc文件，转换为匹配用的像素坐标格式
% input:
% path  .vc文件存放的路径及文件名（含后缀）
% output:
% uvc    1*frameN       cell型，单个镜头所有帧的像素坐标信息
% uvc{}  1*2markerN     每个元素为像素坐标uv信息，
%                       按[ u1 v1 u2 v2 u3 v3 ...]顺序排列，表示第1、2、3个点的u坐标和v坐标

fid=fopen(path,'r')  ;
% fid=fopen(['C:\Users\Boat\Desktop\1.0 vcbuild\VCFiles\walk\walk.vc',int2str(icam+1)],'r')  ;%icam+1第一个镜头标定不好，从2号镜头开始
[si,num] = fread(fid,[1,inf],'uint8');
fclose(fid);

% flag = 0 ;%debug
i=1 ;
while 1
    if si(i)==255 && si(i+1)==2
        nframe = si(i+2) + si(i+3)*16^2 ; %帧号
        i = i+4 ;
        pn = si(i) ; %point_number 这一帧的点数,不超过255个点
        i = i+4 ;
        if pn<1, uvc{nframe}=[NaN, NaN] ; end
        for k = 1:pn
            uvc{nframe}(k*2-1) = readvc_bin2dec([si(i+3), si(i+2), si(i+1), si(i)]) ;
            i = i+4 ;
            uvc{nframe}(k*2) = readvc_bin2dec([si(i+3), si(i+2), si(i+1), si(i)]) ;
            i = i+8 ;
        end
%         flag = flag+1 ;%debug
    end
    i = i+1 ;
     if i+1>length(si), break; end
%     if flag>3, break; end %debug
end

end

function dec = readvc_bin2dec(si)

% 将镜头数据.vc文件转换为匹配可用的像素坐标时用
% 将4个连续的10进制整数转换为十六进制浮点型再变回一个十进制的实数
% si 1*4 10进制整数

si = dec2hex(si) ;
size2 = size(si,2) ;
if size2 ~= 2
    if size2 == 1
        si = [int2str(zeros(4,1)), si] ;
    else
        error ;
    end
end
si = si' ;
si = si(:) ;
si = si' ;
dec = typecast(uint32(hex2dec(si)),'single') ; %point_number 这一帧的点数

end

