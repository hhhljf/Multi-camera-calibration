function inexI = movedL(xf, Luv, L)

% 将整个坐标系平移旋转到以L框架为标准的坐标系上
% input：
% xf    10camN*1    LM优化后的所有镜头的内外参
% Luv   4*2camN     L型框架的4个点在所有镜头中的像素坐标，每4*2个矩阵为一个镜头的uv像素点坐标
%                   像素坐标uv的排列顺序可以不与xyz中的点严格对应。 
% L     4*3         L型框架在世界坐标系下的坐标，通常为
%                   L=[   0	  0	0
%                       200	  0	0
%                       600	  0	0
%                         0	400	0 ];
% 整个过程分为两步
% 1.平移，按当前内外参数重建出L型框架的坐标与实际设定的L型框架坐标比较，
%         将重建出的L型框架坐标原点平移到设定的原点，所有镜头的位移参数随之平移即可；
% 2.旋转，由于当前各镜头已旋转了一定角度，直接在原空间中旋转镜头不易将旋转过程描述为各镜头的欧拉角旋转，
%         重新计算旋转参数更为简单直接，旋转过程描述为：先随L型框架旋转exRd到设定位置，再按各自的旋转参数旋转exR0，
%         整合两次旋转过程exR=exR0*exRd（此处exR0和exRd的位置不可交换），以整体旋转矩阵exR计算旋转参数.

% Luv = [1027	698	908.0999756	634	1023.400024	424	942.5	499.5	850.4000244	609	760	505	935.125	506.5	949	428.5	1118	587.5	852.5	498.5	932.555542	484	1007	528	995.375	495.5	1028.5	516	854.333313	662.5	955.833313	578.5
% 815.75	708.5	825.5	671	924.5	434.5	1046.214233	506	810.666687	631.5	845.5	538	959.5	542.5	935.7999878	470	945.416687	606.5	1110.666626	519.5	887.611084	556	916.5	547	969.7000122	510	960.1428833	527	991.083313	681.5	1003.75	605.5
% 943.5	730	870.7999878	690	935.625	452.5	934.083313	518.5	942.416687	634.5	891.2999878	555	972.5999756	562	811.7999878	489	1059.916626	614.5	972.833313	550.5	962.055542	570	941.4000244	564	1085.5	521.5	818.875	548.5	1062.099976	691	1029.5	619.5
% 898.666687	747.5	967.2000122	729	960.2000122	492	916.166687	561.5	725.583313	678.5	805.833313	595.5	863.5	575	928.416687	492.5	1028.5	629	900.833313	566.5	1114.400024	599.5	994.9000244	601	914.7999878	540	1081	590.5	1021.928589	767	942.666687	636.5
% ];
% L = [0 0 0;
%     200 0 0 ;
%     600 0 0 ;
%     0 400 0 ;
%     ] ;

% load input\xf
camN = length(xf)/10 ;
inexI = reshape(xf,10,camN) ;  % 将 1*10n 的向量变为 10*n 的矩阵
inexI = inexI' ;

M = [];
for i = 1:camN
    temp = buildM(inexI(i,:)) ;
    M = [M temp'];
end
M0=M' ;

% 将Luv按照 0 1 2 3 号点的顺序排序
for i = 1:camN
    Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
end

% 重建出的L型框架的坐标
pwc = zeros(4,3) ;
for k = 1:4
    pwc(k,:) = rebulid_3D_UnspecCam_LM2(Luv(k,:),M0)  ;
end
tran = pwc(1,:) ; %偏移距离
pwc = pwc - repmat(pwc(1,:),4,1) ; 

% 用反对称矩阵和罗德里格矩阵求解旋转矩阵
A = zeros(3*4,3) ; 
for i = 1:4
    A(3*i-2:3*i,:) = [  0   L(i,3)+pwc(i,3)   L(i,2)+pwc(i,2) 
                     L(i,3)+pwc(i,3)      0   -L(i,1)-pwc(i,1)
                     L(i,2)+pwc(i,2)  L(i,1)+pwc(i,1)     0 ] ;
    B(3*i-2:3*i) = [L(i,1)-pwc(i,1); L(i,2)-pwc(i,2); L(i,3)+pwc(i,3)] ;
end
B=B' ;
temp = A\B ;
S = [0 -temp(3) -temp(2); temp(3) 0 -temp(1); temp(2) temp(1) 0];
exRd = (eye(3,3)+S)*inv((eye(3,3)-S)) ;  % 旋转矩阵 R=(I+S)*(I-S)^(-1)

% 计算位移参数
inexI(:,5:7) = inexI(:,5:7) - repmat(tran,camN,1) ;
temp = exRd'*inexI(:,5:7)' ;
inexI(:,5:7) = temp' ;

% 计算旋转参数
temp = inexI(:,8:10) ;
temp = temp*pi/180 ;
for i = 1:size(temp,1)
    rx = temp(i,1); ry = temp(i,2); rz = temp(i,3); 
    exR0=[cos(ry)*cos(rz) cos(ry)*sin(rz) -sin(ry)
    sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz) sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz) sin(rx)*cos(ry)
    cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz) cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz) cos(rx)*cos(ry)
    ];
    temp(i,:) = exR2rxryrz(exR0*exRd) ;  %这里 exR0*exRd 的顺序不能变，按先随L型框架旋转exRd，再按各自的旋转参数旋转exR0的顺序
end
inexI(:,8:10) = temp*180/pi ;

end

function r = exR2rxryrz(exR)
% 将旋转矩阵exR按旋转顺序 rz -> ry -> rx，分解为rx,ry,rz
% input:
% exR   3*3     旋转矩阵，顺序 rz -> ry -> rx
% output:
% r     1*3     [rx ry rz]

ry = -asin(exR(1,3)) ;
rx = acos(exR(3,3)/cos(ry));
% if ( rx_hat>0 ) rx_hat=-rx_hat ; end  % 镜头都是俯视，所以只会绕x往负方向旋转. 

if sign(sin(rx)*cos(ry)) ~= sign(exR(2,3))
    ry = sign(ry)*pi - ry ; 
    rx = acos(exR(3,3)/cos(ry));
%     if ( rx_hat>0 ) rx_hat=-rx_hat ; end 
end

rz = asin(exR(1,2)/cos(ry));

if sign(cos(ry)*cos(rz)) ~= sign(exR(1,1))
    rz = sign(rz)*pi - rz ;
end

r = [rx ry rz ] ;

end

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

