function orderuv = Lframe(uv) 
%�������ܣ�����Lframe�ϵ�4���㣬�����������г�Ϊ 0 1 2 3 �ŵ�˳�� .uv �� orderuv ����4*2��ʽ
% 0 1 2 ��x���ϣ�3 ��y���ϡ�0Ϊ����ԭ�㡣

% ����˼�룺������ȷ��һ��ֱ�ߡ�����������������ֱ�߾���ܽ�����ֵ����Ϊ��5�����Ե��ڣ���Ϊ��һ��ֱ���ϵġ�����ֱ�ߵ��ϵ���3�ŵ㡣
% ���� 0 1 2 ��֮��ľ��롣��������� 0 2 �㡣�ɴ�ȷ��������һ������1�ŵ㡣 ����1�ŵ���������0�ŵ㣬����Զ����2�ŵ㡣

% Ϊ�˱�����д�����е��ȷ����
uv0=uv(1,:) ; 
uv1=uv(2,:) ;
uv2=uv(3,:) ;
uv3=uv(4,:) ;

orderuv=zeros(4,2) ;
% d ��ʱ��������ʾ����

threshold =5 ;% ��ֵ�趨Ϊ5

%�� 0 1 ����ֱ�ߣ������ж�

d1=disPtoLine(uv0,[uv1 uv2]);
d2=disPtoLine(uv0,[uv1 uv3]);
d3=disPtoLine(uv0,[uv2 uv3]);

if d1>threshold && d2>threshold && d3>threshold
    orderuv(4,:)=uv0;
    orderuv(1:3,:)=Marker_sort([uv3;uv1;uv2]) ;
%     ��������
end


if d1<threshold
    orderuv(4,:)=uv3;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv2]) ;
%     ��������
end

    
    
if d2<threshold
    orderuv(4,:)=uv2;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv3]) ;
%     ��������
    
end    
    
if d3<threshold
    orderuv(4,:)=uv1;
    orderuv(1:3,:)=Marker_sort([uv0;uv3;uv2]) ;
%     ��������
    
end

end %Lframe end



function New_marker=Marker_sort(marker)
%1,2,3�ŵ�����
New_marker = zeros(3,2);
d1_2=sqrt((marker(2,1)-marker(1,1))^2+(marker(2,2)-marker(1,2))^2); %1,2�ŵ�֮��ľ���
d2_3=sqrt((marker(2,1)-marker(3,1))^2+(marker(2,2)-marker(3,2))^2); %2,3�ŵ�֮��ľ���
d1_3=sqrt((marker(1,1)-marker(3,1))^2+(marker(1,2)-marker(3,2))^2); %1,3�ŵ�֮��ľ���
d = [ d2_3 d1_3 d1_2 ];
index = find(d==max(d));
New_marker(2,:) = marker( index,:) ; %��������� 0 2 �㡣�ɴ�ȷ��������һ������1�ŵ㡣
index = find(d<max(d)) ;
New_marker(1,:) = marker(find(d==max(d(index))),:) ; %����1�ŵ���������0�ŵ㣬
New_marker(3,:) = marker(find(d==min(d(index))),:) ; %����Զ����2�ŵ㡣
end %Marker_sort end

function [ d ] = disPtoLine( P, CorPts ) 
%�������ܣ�����㵽ֱ�ߵľ��루���������ά�ģ�
% ���룺corpts ��u1 v1 u2 v2����ʽ������ȷ��һ��ֱ�ߡ�
% P ��u v����ʽ
% �����d ��ʾ�㵽ֱ�ߵľ���

% Detailed explanation goes here 
if length(P) == 2
    l = [ CorPts(1) - CorPts(3), CorPts(2) - CorPts(4), 0 ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), 0 ]; 
elseif length(P) == 3
    l = [ CorPts(1) - CorPts(4), CorPts(2) - CorPts(5), CorPts(3) - CorPts(6) ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), P(3) - CorPts(3) ]; 
else
    error('������������~') ;
end
tem = cross(pl, l); 
d = norm( tem ) / norm( l ); 

end %disPtoLine end


