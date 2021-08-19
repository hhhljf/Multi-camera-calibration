function orderuv = Lframe(uv) 
%函数功能：输入Lframe上的4个点，进行排序。排列成为 0 1 2 3 号的顺序 .uv 及 orderuv 都是4*2格式
% 0 1 2 在x轴上；3 在y轴上。0为坐标原点。

% 程序思想：两个点确定一条直线。其他点如果距离这个直线距离很近（阈值，认为是5，可以调节）认为是一条直线上的。不在直线的上点是3号点。
% 计算 0 1 2 点之间的距离。距离最长的是 0 2 点。由此确定，另外一个点是1号点。 距离1号点距离近的是0号点，距离远的是2号点。

% 为了便于书写，进行点的确定。
uv0=uv(1,:) ; 
uv1=uv(2,:) ;
uv2=uv(3,:) ;
uv3=uv(4,:) ;

orderuv=zeros(4,2) ;
% d 临时变量，表示距离

threshold =5 ;% 阈值设定为5

%对 0 1 建立直线，进行判定

d1=disPtoLine(uv0,[uv1 uv2]);
d2=disPtoLine(uv0,[uv1 uv3]);
d3=disPtoLine(uv0,[uv2 uv3]);

if d1>threshold && d2>threshold && d3>threshold
    orderuv(4,:)=uv0;
    orderuv(1:3,:)=Marker_sort([uv3;uv1;uv2]) ;
%     进行排序
end


if d1<threshold
    orderuv(4,:)=uv3;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv2]) ;
%     进行排序
end

    
    
if d2<threshold
    orderuv(4,:)=uv2;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv3]) ;
%     进行排序
    
end    
    
if d3<threshold
    orderuv(4,:)=uv1;
    orderuv(1:3,:)=Marker_sort([uv0;uv3;uv2]) ;
%     进行排序
    
end

end %Lframe end



function New_marker=Marker_sort(marker)
%1,2,3号点排序
New_marker = zeros(3,2);
d1_2=sqrt((marker(2,1)-marker(1,1))^2+(marker(2,2)-marker(1,2))^2); %1,2号点之间的距离
d2_3=sqrt((marker(2,1)-marker(3,1))^2+(marker(2,2)-marker(3,2))^2); %2,3号点之间的距离
d1_3=sqrt((marker(1,1)-marker(3,1))^2+(marker(1,2)-marker(3,2))^2); %1,3号点之间的距离
d = [ d2_3 d1_3 d1_2 ];
index = find(d==max(d));
New_marker(2,:) = marker( index,:) ; %距离最长的是 0 2 点。由此确定，另外一个点是1号点。
index = find(d<max(d)) ;
New_marker(1,:) = marker(find(d==max(d(index))),:) ; %距离1号点距离近的是0号点，
New_marker(3,:) = marker(find(d==min(d(index))),:) ; %距离远的是2号点。
end %Marker_sort end

function [ d ] = disPtoLine( P, CorPts ) 
%函数功能：计算点到直线的距离（下面的是三维的）
% 输入：corpts （u1 v1 u2 v2）格式，用来确定一条直线。
% P （u v）格式
% 输出：d 表示点到直线的距离

% Detailed explanation goes here 
if length(P) == 2
    l = [ CorPts(1) - CorPts(3), CorPts(2) - CorPts(4), 0 ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), 0 ]; 
elseif length(P) == 3
    l = [ CorPts(1) - CorPts(4), CorPts(2) - CorPts(5), CorPts(3) - CorPts(6) ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), P(3) - CorPts(3) ]; 
else
    error('坐标数量错误~') ;
end
tem = cross(pl, l); 
d = norm( tem ) / norm( l ); 

end %disPtoLine end


