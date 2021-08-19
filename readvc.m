function uvc = readvc(path)

% ��ȡ����.vc�ļ���ת��Ϊƥ���õ����������ʽ
% input:
% path  .vc�ļ���ŵ�·�����ļ���������׺��
% output:
% uvc    1*frameN       cell�ͣ�������ͷ����֡������������Ϣ
% uvc{}  1*2markerN     ÿ��Ԫ��Ϊ��������uv��Ϣ��
%                       ��[ u1 v1 u2 v2 u3 v3 ...]˳�����У���ʾ��1��2��3�����u�����v����

fid=fopen(path,'r')  ;
% fid=fopen(['C:\Users\Boat\Desktop\1.0 vcbuild\VCFiles\walk\walk.vc',int2str(icam+1)],'r')  ;%icam+1��һ����ͷ�궨���ã���2�ž�ͷ��ʼ
[si,num] = fread(fid,[1,inf],'uint8');
fclose(fid);

% flag = 0 ;%debug
i=1 ;
while 1
    if si(i)==255 && si(i+1)==2
        nframe = si(i+2) + si(i+3)*16^2 ; %֡��
        i = i+4 ;
        pn = si(i) ; %point_number ��һ֡�ĵ���,������255����
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

% ����ͷ����.vc�ļ�ת��Ϊƥ����õ���������ʱ��
% ��4��������10��������ת��Ϊʮ�����Ƹ������ٱ��һ��ʮ���Ƶ�ʵ��
% si 1*4 10��������

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
dec = typecast(uint32(hex2dec(si)),'single') ; %point_number ��һ֡�ĵ���

end

