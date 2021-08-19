
function exI = LFrameCalibration_03(inI,uv)
% This design is the property of NOKOV.  Publication of this
% design is not authorized without written consent from NOKOV.
% 
% Disclaimer:
%    NOKOV makes no warranty for the use of this code or design.
%----------------------------------------------------------------
%
% Create Date:         May 09, 2016
% Design Name:         NOKOV 1.1.0.1
% Author:               Boat
% Tool versions:      MATLAB 7.14.0.739 (R2012a)
% Revision:            Aug 11, 2016: 1.01 Initial version
%----------------------------------------------------------------
%
% Function Description:    
% 在已知单镜头内参和L型框架信息的情形下，计算单镜头外参初值。
% 需要调用外部函数 LMFnlsq2.m,程序运行中会生成文件input\xyzC.txt .
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input Arguments:
% inI 1*4 内参向量 inI = [f dx u0 v0]
%         f  镜头焦距
%         dx 像素点在物理成像面x轴上的尺寸比例
%         u0 成像面坐标系的主点在像素坐标系下的横坐标
%         v0 成像面坐标系的主点在像素坐标系下的纵坐标
%         其他的内参 dy固定为0.0048 单位mm; ac固定为pi/2。
% uv  4*2 L型框架的4个点在单镜头中的像素坐标
%         L型框架的世界坐标固定为
%         xyz=[0	0	0
%              200	0	0
%              600	0	0
%              0	400	0 ];
%         像素坐标uv的排列顺序可以不与xyz中的点严格对应。 
%
% Output Arguments:
% exI 1*6 外参向量 exI = [tx ty tz rx ry rz]
%                  tx  镜头在世界坐标系下的x轴坐标
%                  ty  镜头在世界坐标系下的y轴坐标
%                  tz  镜头在世界坐标系下的z轴坐标
%                  rx  镜头从初始位置（世界坐标原点）开始，绕x轴旋转的欧拉角角度
%                  ry  镜头从初始位置（世界坐标原点）开始，绕y轴旋转的欧拉角角度
%                  rz  镜头从初始位置（世界坐标原点）开始，绕z轴旋转的欧拉角角度
%                  欧拉角旋转顺序 rz -> ry -> rx
%
%% Example :
%
inI = [8.628	0.0048	980	500] ;
uv  = [ 770.1, 607.12
       859.14, 617.79
       1036.5, 639.05
       791.65, 534.51];
%    
% exI = LFrameCalibration_03
%
% exI =
% 
%        327.84        -3597       1831.7      -112.98       -7.005       1.2846
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 历史版本
% 2016.03.25 waican_L4pt_v1.m 
%            在waican_L4pt.m上的改动：用了新的方法计算P4点在相机坐标系下的位置
%
% 2016.04.05 waican_L4pt_v2.m 
%            改变了计算旋转角度的判断方式。
%            角度错误主要体现在三角函数的符号上，因此只要判断用rx,ry,rz计算出的矩阵各元素符号是否与exR_hat相同即可。
%            将 Lframe.m 和 disPtoLine.m 两个函数的内容放到了这个函数中，不需要再调用外部 .m 文件。
%
% 2016.04.11 waican_L4pt_v3.m
%            1)在得到4个点在相机坐标系下的坐标后，用L型框架的约束条件做一次LM优化，调整4个点的坐标，目标函数为lm4pt.m .
%              在得到成像面坐标xyzC后，将其保存在input\xyzC.txt中，供lm4pt.m调用。
%            2)因为位移参数的算法对像素坐标uv的敏感度太高，而旋转参数所受影响较小。
%              所以在得到旋转参数rx,ry,rz后，用旋转参数反求位移参数tx,ty,tz，减小像素坐标的偏差对位移参数的影响。
%
% 2016.05.06 LFrameCalibration.m
%            将需要外部调用的lm4pt.m整合到了此函数下，删除了部分debug注释；
%            整个标定过程整理为5部分：
%            1.输入数据的预处理
%            2.求解4个点在相机坐标系下的坐标
%              2.1 计算共线的3个点的相机坐标
%              2.2 计算第4号点的坐标
%              2.3 以L型框架的约束条件对4个点的坐标进行优化
%            3.将世界坐标的原点平移到相机坐标原点位置
%            4.计算旋转参数
%            5.反求位移参数
%
% 2016.08.11 LFrameCalibration_01.m
%            2.2步，计算4号点在相机坐标系下的坐标中，
%            若光心位于过Pw4且与Pw2Pw1垂直的平面上，会出现Pw4_1和Pw4_2垂直度相同（均为90度）的情况，
%            因此，当光心在过Pw4且与Pw2Pw1垂直的平面附近时，只用垂直度判断容易出错，需要添加判断条件：
%            若2D平面（像素平面或归一化平面）中P4的位置在P1之下，则取Pw4_1和Pw4_2中离光心O较近的点；
%            反之，若D平面中P4的位置在P1之上，则取Pw4_1和Pw4_2中离光心O较近的点。
%            这个条件是基于摄像头总是在世界坐标系的Z=0平面之上、俯视整个拍摄区域的，其他情形需要调整条件。
% 
% 2016.10.20 LFrameCalibration_02.m
%            利用LM算法优化4个点在相机坐标系下的坐标部分，将LMFnlsq2修改为L型标定专用的LM算法函数Calibration_Lframe_LM2，
%            直接使用归一化平面上的投影坐标，不再保存在外部文件中，并对应修改了lm4pt的输入参数。
% 
% 2016.10.22 LFrameCalibration_03.m
%            优化2.1计算共线的3个点的相机坐标，2）1,2,3号点共线，且有比例关系 3*P1P2 = P1P3
%            整合约束条件，将系数矩阵A的行数减少为9行。
%            提高适用性，允许自由设定L型框架上四个点的世界坐标。
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
% format long

% ---------------------------------------------------------------------
% ----------------------- 1.输入数据的预处理 ---------------------------

[uv_l,uv_c]=size(uv); %求取uv的行数、列数
[inI_l,inI_c]=size(inI);   %求取inI的行数、列数
 if inI_l~=1 || inI_c~=4
     error('输入内参inI的格式为%d行%d列，应为1行4列\n',inI_l,inI_c);
 end 
 if uv_l~=4 || uv_c~=2
     error('输入L型框架像素坐标uv的格式为%d行%d列，应为4行2列\n',inI_l,inI_c);
 end 

xyz=[0	0	0 
200	0	0
600	0	0
0	400	0
]; % 规定好的L型框架世界坐标
uv = Lframe(uv) ; %对uv进行排序，使之与xyz坐标匹配

ac=pi/2;
f=inI(1);
dx=inI(2);
dy=0.0048;
fx=f/dx;
fy=f/dy;
u0=inI(3);
v0=inI(4);

in=[fx -fx/tan(ac) u0
    0 fy/sin(ac) v0
    0 0 1]; 

UV = [uv,ones(size(uv,1),1)];
xyzC = (in\(UV'))' ;  % 将像素上的点转换为相机坐标系归一化平面上的点
% dlmwrite('input\xyzC.txt',xyzC); %供函数lm4pt在LM优化时调用

% --------------------- 1.输入数据的预处理 end -------------------------
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% ---- 2.利用小孔成像原理构建线性方程组，求解4个点在相机坐标系下的坐标 ----

% 2.1 计算共线的3个点的相机坐标
% 约束条件：
% 1） 光心、相机投影平面上的点、世界坐标上对应的点 在同一直线上
% 2） 1,2,3号点共线，且有比例关系 3*P1P2 = P1P3
tic
L13 = sqrt( (xyz(3,1)-xyz(1,1))^2 + (xyz(3,2)-xyz(1,2))^2 + (xyz(3,3)-xyz(1,3))^2 ) ; %世界坐标系下1号点与3号点间的距离
L12 = sqrt( (xyz(2,1)-xyz(1,1))^2 + (xyz(2,2)-xyz(1,2))^2 + (xyz(2,3)-xyz(1,3))^2 ) ; %世界坐标系下1号点与2号点间的距离
A = zeros(9,9) ;
for i=1:3
    A(i*2-1, 1+3*(i-1)) = 1 ;
    A(i*2-1, 3*i)       = -xyzC(i,1) ;
    A(i*2  , 2+3*(i-1)) = 1 ;
    A(i*2  , 3*i)       = -xyzC(i,2) ;
    A(6+i  , i)         = L13 - L12 ;
    A(6+i  , i+3)       = -L13 ;
    A(6+i  , i+6)       = L12 ;
end
% rref(A) ; % debug
[U S V] = svd(A) ;
pw1=V(1:3,9)';  
pw2=V(4:6,9)';  
pw3=V(7:9,9)';  
bili =  L12 / sqrt((pw2-pw1)*(pw2-pw1)') ;  % 用P1P2之间的距离200，计算 解出的坐标 与 实际坐标 之间的比例
pw1 = pw1 * bili ;  % P1点在相机坐标系中的坐标
pw2 = pw2 * bili ;  % P2点在相机坐标系中的坐标
pw3 = pw3 * bili ;  % P3点在相机坐标系中的坐标

% pw1(3)<0时，求到了相机背后的点，所有坐标反向即可.
if pw1(3) < 0
    pw1 = -pw1 ;
    pw2 = -pw2 ;
    pw3 = -pw3 ;
end
exT_hat = pw1 ; % pw1 的坐标就是 exT 的坐标

% 2.2 计算第4号点的坐标
% pw4 坐标的计算方法 
% 1.找出射线OP4上到P1的距离为400的点，正常情况下有两个.（若不存在，则P4为P1到OP4的垂点.）
% 2.过这两个点与P1做两条直线，与P1P2垂直度较高的线对应的点作为P4点的位置.
% 3.若光心在过P4与P1P2垂直的平面附近，另外判断

L14 = sqrt( (xyz(4,1)-xyz(1,1))^2 + (xyz(4,2)-xyz(1,2))^2 + (xyz(4,3)-xyz(1,3))^2 ) ; %世界坐标系下1号点与4号点间的距离
d14 = disPtoLine( pw1, [0 0 0 xyzC(4,:)] ) ; % P1到直线OP4的距离
theta4x = acos(xyzC(4,1)/norm(xyzC(4,:))) ; % 直线OP4与x轴的夹角
theta4y = acos(xyzC(4,2)/norm(xyzC(4,:))) ; % 直线OP4与y轴的夹角
theta4z = acos(xyzC(4,3)/norm(xyzC(4,:))) ; % 直线OP4与z轴的夹角

theta = acos( pw1*xyzC(4,:)'/(norm(pw1)*norm(xyzC(4,:))) ) ;   % OP1与OP4的夹角
l4_1 = cos(theta)*norm(pw1) - sqrt(max(0,L14^2-d14^2)) ; %P4在直线OP4上有可能的位置1
l4_2 = cos(theta)*norm(pw1) + sqrt(max(0,L14^2-d14^2)) ; %P4在直线OP4上有可能的位置2
pw4_1 = l4_1*cos([theta4x theta4y theta4z]) ; %P4位置1的点的坐标
pw4_2 = l4_2*cos([theta4x theta4y theta4z]) ; %P4位置2的点的坐标
theta4_1 = acos( abs((pw1-pw4_1)*(pw1-pw2)')/(norm(pw1-pw4_1)*norm(pw1-pw2)) ) ;  %P1P4位置1的点与P1P2的夹角
theta4_2 = acos( abs((pw1-pw4_2)*(pw1-pw2)')/(norm(pw1-pw4_2)*norm(pw1-pw2)) ) ;  %P1P4位置2的点与P1P2的夹角
if abs(pi/2-abs(theta4_1)) < abs(pi/2-abs(theta4_2)) % 与P1P2垂直度较高的作为P4点的位置
    pw4 = pw4_1 ;
else
    pw4 = pw4_2 ;
end
toc
% 若光心在过P4与P1P2垂直的平面附近，会造成P4可能位置的两个点pw4_1,pw4_2垂直度均较高的情况，
% 此时若只用垂直度判断容易出错，需要添加条件.
% 若2D平面中P4的像在P1之下，则取pw4_1和pw4_2中Z坐标较小的那一个；反之，
% 若2D平面中P4的像在P1之上，则取pw4_1和pw4_2中Z坐标较大的那一个.
% ！！！！注意，此处假设摄像头总是在世界坐标系的Z=0平面之上、俯视整个拍摄区域的！！！！
if abs(theta4_1*180/pi)>80 & abs(theta4_2*180/pi)>80   %此处80为判断垂直的阈值，可以修改
    if uv(4,2) > uv(1,2)
        if pw4_1(3)<pw4_2(3)
            pw4 = pw4_1;
        else
            pw4 = pw4_2 ;
        end
    else
        if pw4_1(3)<pw4_2(3)
            pw4 = pw4_2;
        else
            pw4 = pw4_1 ;
        end
    end
end

% % 若摄像头在世界坐标系的Z=0平面之下、仰视整个拍摄区域
% if abs(theta4_1*180/pi)>80 & abs(theta4_2*180/pi)>80   %此处80为判断垂直的阈值，可以修改
%     if uv(4,2) > uv(1,2)
%         if pw4_1(3)<pw4_2(3),  pw4 = pw4_2;
%         else pw4 = pw4_1 ;
%         end
%     else
%         if pw4_1(3)<pw4_2(3),  pw4 = pw4_1;
%         else pw4 = pw4_2 ;
%         end
%     end
% end


% 2.3 以L型框架的约束条件对4个点的坐标进行优化
pw0 = [pw1,pw2,pw3,pw4];

xf = Calibration_Lframe_LM2(pw0, xyzC, xyz); %LM

pw1 = xf(1:3)' ;
pw2 = xf(4:6)' ;
pw3 = xf(7:9)' ;
pw4 = xf(10:12)' ;

% ----- 2.利用小孔成像原理构建方程组，4个点在相机坐标系下的坐标  end  -----
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% ------------------------------ 3.平移 -------------------------------

pw = [pw1;pw2;pw3;pw4]; %四个点在相机坐标系下的坐标
pwc = pw - repmat(pw1,4,1) ; %将世界坐标的原点(pw1)平移到相机坐标原点位置

% ---------------------------- 3.平移 end -----------------------------
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% ----------------------- 4.计算旋转参数rx ry rz -----------------------
 
%  用反对称矩阵和罗德里格矩阵求解旋转矩阵
A = zeros(3*4,3) ; 
for i = 1:4
    A(3*i-2:3*i,:) = [  0   xyz(i,3)+pwc(i,3)   xyz(i,2)+pwc(i,2) 
                     xyz(i,3)+pwc(i,3)      0   -xyz(i,1)-pwc(i,1)
                     xyz(i,2)+pwc(i,2)  xyz(i,1)+pwc(i,1)     0 ] ;
    B(3*i-2:3*i) = [xyz(i,1)-pwc(i,1); xyz(i,2)-pwc(i,2); xyz(i,3)-pwc(i,3)] ;
end
B=B' ;
temp = A\B ;
S = [0 -temp(3) -temp(2); temp(3) 0 -temp(1); temp(2) temp(1) 0];
exR_hat = (eye(3,3)+S)/(eye(3,3)-S) ;  % 旋转矩阵 R=(I+S)*(I-S)^(-1)

% 用反对称矩阵和罗德里格矩阵求解旋转矩阵 end 

% 计算旋转角度
ry_hat = -asin(exR_hat(1,3)) ;
rx_hat = acos(exR_hat(3,3)/cos(ry_hat));
if ( rx_hat>0 ) rx_hat=-rx_hat ; end  % 镜头都是正放，所以只会绕x往负方向旋转. 
                                      % 若镜头倒放， if ( rx_hat<0 ) rx_hat=-rx_hat ; end

if sign(sin(rx_hat)*cos(ry_hat)) ~= sign(exR_hat(2,3))
    ry_hat = sign(ry_hat)*pi - ry_hat ; 
    rx_hat = acos(exR_hat(3,3)/cos(ry_hat));
    if ( rx_hat>0 ) rx_hat=-rx_hat ; end % 镜头都是正放，所以只会绕x往负方向旋转. 
                                         % 若镜头倒放， if ( rx_hat<0 ) rx_hat=-rx_hat ; end
end

rz_hat = asin(exR_hat(1,2)/cos(ry_hat));

if sign(cos(ry_hat)*cos(rz_hat)) ~= sign(exR_hat(1,1))
    rz_hat = sign(rz_hat)*pi - rz_hat ;
end
% 计算旋转角度 end

% --------------------- 4.计算旋转参数rx ry rz end ---------------------
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% -------- 5.用旋转参数反求位移参数，减小像素偏差对位移参数的影响 ---------

rx = rx_hat;
ry = ry_hat;
rz = rz_hat ;
exR0=[cos(ry)*cos(rz) cos(ry)*sin(rz) -sin(ry)
    sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz) sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz) sin(rx)*cos(ry)
    cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz) cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz) cos(rx)*cos(ry)
    ];

A = zeros(8,3); 
b = zeros(8,1);
for i = 1:4
    A(2*i-1:2*i,:) =[ 1 0 -xyzC(i,1) 
                      0 1 -xyzC(i,2)] ;
    b(2*i-1:2*i) = ( exR0(1:2,:) - xyzC(i,1:2)'*exR0(3,:) ) * xyz(i,:)' ;
end
temp = A\b ;
exI = [ (exR0'*temp(1:3))' rx_hat*180/pi ry_hat*180/pi rz_hat*180/pi] ;

% ------ 5.用旋转参数反求位移参数，减小像素偏差对位移参数的影响 end -------
% ---------------------------------------------------------------------


end %LFramCalibration end

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

% function d = lm4pt(x, xyzC)
% % 在得到4个点在相机坐标系下的坐标后，用L型框架的约束条件做一次LM优化的目标函数，用于调整4个点的坐标。
% % 所用到的约束条件有：
% % 1. 4个点之间的距离；
% % 2. 1）四个点都在OP_i射线上；2）1,2,3号点比例共线；
% % 3. P4P1 垂直于 P1P2, P1P3, P2P3 .
% % 一共26个方程，在目标函数中，将第1部分约束权值加大，可以得到比以前好的效果。
% % 输入：x 1*12 依次为1,2,3,4号点的xyz坐标
% % 输出：最小二乘方法中的求和项，理想情况下有 d(i)=0 .
% 
% % xyzC =[
% %    0.067985187966913   0.195717726274982   1.000000000000000
% %    0.066622629957160   0.157654787086245   1.000000000000000
% %    0.064295466900847   0.092643660162677   1.000000000000000
% %   -0.059236541965606   0.180997535670956   1.000000000000000
% % ];
% 
% % xyzC = importdata('input\xyzC.txt') ; 
% 
% x = x(:) ; 
% P = reshape(x,3,4) ;  
% P = P' ; % P 4*3, 每一行代表一个点的xyz坐标
% 
% d = zeros(20,1) ; 
% 
% % 1. 四个点两两之间的距离
% d(1) = sqrt((P(2,:)-P(1,:))*(P(2,:)-P(1,:)).') - 200 ; % 1 2 号点之间的距离为200
% d(2) = sqrt((P(3,:)-P(1,:))*(P(3,:)-P(1,:)).') - 600 ; % 1 3 号点之间的距离为600
% d(3) = sqrt((P(4,:)-P(1,:))*(P(4,:)-P(1,:)).') - 400 ; % 1 4 号点之间的距离为400
% d(4) = sqrt((P(3,:)-P(2,:))*(P(3,:)-P(2,:)).') - 400 ; % 2 3 号点之间的距离为400
% d(5) = sqrt((P(4,:)-P(2,:))*(P(4,:)-P(2,:)).') - sqrt(200^2+400^2) ; % 2 4 号点之间的距离为sqrt(200^2+400^2)
% d(6) = sqrt((P(4,:)-P(3,:))*(P(4,:)-P(3,:)).') - sqrt(600^2+400^2) ; % 3 4 号点之间的距离为sqrt(600^2+400^2)
% 
% % 2.1 四个点都在OP_i射线上
% A = zeros(11,12) ;
% 
% for i=1:4
%     A(i*2-1, 1+3*(i-1)) = 1 ;
%     A(i*2-1, 3*i)       = -xyzC(i,1) ;
%     A(i*2  , 2+3*(i-1)) = 1 ;
%     A(i*2  , 3*i)       = -xyzC(i,2) ;
% end
% % 2.2 1,2,3号点比例共线
% A(9:11,1:3) = 2*eye(3) ;
% A(9:11,4:6) = -3*eye(3) ;
% A(9:11,7:9) = eye(3) ;
% 
% d(7:17) = A*x ;
% 
% % 3. P4P1 垂直于 P1P2, P1P3, P2P3
% d(18) = ((P(4,:)-P(1,:))*(P(2,:)-P(1,:))')/(norm(P(4,:)-P(1,:))*norm(P(2,:)-P(1,:))) ;
% d(19) = ((P(4,:)-P(1,:))*(P(3,:)-P(1,:))')/(norm(P(4,:)-P(1,:))*norm(P(3,:)-P(1,:))) ;
% d(20) = ((P(4,:)-P(1,:))*(P(3,:)-P(2,:))')/(norm(P(4,:)-P(1,:))*norm(P(3,:)-P(2,:))) ;
% 
% end % lm4pt end
