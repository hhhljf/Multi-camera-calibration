function M = buildM(inexI)


% 函数功能：输入inexI内外参数的物理数值，输出M矩阵，一共1行10列。
% % % % % 构建M矩阵，进行归一化操作% % % % % 

% % % 利用内外参数的实际数值，构建内外参矩阵和M矩阵
% inexI,表示内参输入矩阵，1行4列，和 表示外参输入矩阵，1行6列。一共一行10列。
% 内参f dx u0 v0  处于inexI的1-4位置。 为使用方便，dy 进行固定，0.0048 单位mm； ac固定pi/2。
% tx,ty,tz,rx,ry,rz,分别表示三个方向上的旋转和平移量。 处于inexI的5-10位置

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% fx,fy,ac,u0 v0,表示内参的五个参数。fx，fy表示焦距，像素表示。ac垂直度，弧度表示
% in 表示内参矩阵，3*3
% ex 表示外参矩阵，3*4 [exR exT]
% exR 外参旋转参数矩阵
% exT 外参平移参数矩阵
% M 表示M矩阵,横排格式1*11
% FM 表示完整M矩阵 3*4
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % %此处输入外参是弧度，需要进行一下转换

% % % % %操作仅仅为便于公式编写及读懂 % % % % % 

ac=pi/2;
f=inexI(1);
dx=inexI(2);
% dx=0.0048;
dy=0.0048;

fx=f/dx;
fy=f/dy;

u0=inexI(3);
v0=inexI(4);

tx=inexI(5);
ty=inexI(6);
tz=inexI(7);
rx=inexI(8)*pi/180; %此处输入外参是弧度，进行一下转换
ry=inexI(9)*pi/180 ;%此处输入外参是弧度，进行一下转换
rz=inexI(10)*pi/180 ;%此处输入外参是弧度，进行一下转换
% % % % %操作仅仅为便于公式编写及读懂 end % % % %

% % % % %构建内参矩阵。公式见说明 % % % % %
% % 此处内参矩阵，可能还可以直接利用参数，后续有可能需要改进
in=[fx -fx/tan(ac) u0
    0 fy/sin(ac) v0  
    0 0 1]; 

% % % % %构建外参矩阵。公式见说明 % % % % % 



% % % XYZ欧拉顺序 Rz * Rz *Rx 
% exR=[cos(rz)*cos(ry) sin(rz)*cos(rx)+cos(rz)*sin(ry)*sin(rx) sin(rz)*sin(rx)-cos(rz)*sin(ry)*cos(rx) 
%     -sin(rz)*cos(ry) cos(rz)*cos(rx)-sin(rz)*sin(ry)*sin(rx) cos(rz)*sin(rx)+sin(rz)*sin(ry)*cos(rx)
%     sin(ry) -cos(ry)*sin(rx) cos(ry)*cos(rx)];
% %     ex公式已经推导验证，正交矩阵性能已经推导测试验证.

% % % ZYX 欧拉顺序....Rx* Rz *Rz 
exR=[cos(ry)*cos(rz) cos(ry)*sin(rz) -sin(ry)
    sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz) sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz) sin(rx)*cos(ry)
    cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz) cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz) cos(rx)*cos(ry)
    ]; % ex公式已经推导验证，正交矩阵性能已经推导测试验证.  
% 
%   




TemexT=[tx;ty;tz];
exT=(-1)*exR*TemexT;

ex=[exR exT];   
 
FM=in*ex;
ss=FM/FM(3,4);
M=[ss(1,:) ss(2,:) ss(3,1) ss(3,2) ss(3,3)];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % 正交矩阵测试程序% % % % % % % 
% 

% r1=exR(:,1)
% r2=exR(:,2)
% r3=exR(:,3)
% r=[r1'*r1 r2'*r2 r3'*r3
%     r1'*r2 r2'*r3 r3'*r1]
% % % % % % % 正交矩阵测试程序end% % % % % % % 


end





