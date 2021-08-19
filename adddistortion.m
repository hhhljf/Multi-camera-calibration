function uvr = adddistortion(uvi, inexI, k)

% 单镜头，在理想的像素坐标uvi上添加畸变，得到实际拍到的像素坐标uvr
% input：
% uvi   1*2n    uv_ideal，畸变校正后n个点理想的像素坐标，每两列为一个点的[u,v]坐标，[u1 v1 u2 v2 ...]
% inexI 1*10    单个镜头的内外参数，格式为[ f dx u0 v0 tx ty tz rx ry rz ]
% k     1*c     c个畸变参数
% output:
% uvr   1*2n    uv_real，实际采集到的n个点的像素坐标，每两列为一个点的[u,v]坐标，[u1 v1 u2 v2 ...]

np = length(uvi)/2 ; %n_point，点的个数
if length(k) ~= 2; error; end % 畸变系数的个数暂定为两个.
u0 = inexI(3) ;
v0 = inexI(4) ;

uvr = zeros(1,2*np) ; % 加上畸变后的像素坐标
for  i = 1:np
    xyi = uvtoim(inexI,uvi(i*2-1:i*2)) ;
    u = uvi(i*2-1) ;
    v = uvi(i*2) ;
%     r2 = norm(xyi) ;
    r2 = xyi(1)^2+xyi(2)^2 ;
    uvr(i*2-1:i*2) = ( [u-u0,0; 0,v-v0]*[1,1;1,1]*[r2,0; 0,r2^2]*k(:) + [u;v] )' ;
end

% % ---parfor---
% uvi = uvi';
% uvi = reshape(uvi, 2, np) ;
% uvi = uvi' ;
% uvr = zeros(np,2) ;
% 
% k=k(:) ;
% parfor i = 1:np
%     xyi = uvtoim(inexI,uvi(i*2-1:i*2)) ;
%     u = uvi(i*2-1) ;
%     v = uvi(i*2) ;
%     r2 = norm(xyi) ;
%     uvr(i,:) = ( [u-u0,0; 0,v-v0]*[1,1;1,1]*[r2,0; 0,r2^2]*k + [u;v] )' ;
% end
% uvr = uvr';
% uvr = uvr(:);
% uvr = uvr' ;

end


function xyi = uvtoim(inexI,uv)

% 函数功能：输入uv，内参数物理数值，输出图像物理坐标xyi矩阵，
% uv,xyi 1*2

% inI,表示内参输入矩阵，1行4列，
% 内参f dx u0 v0  处于inexI的1-4位置。 为使用方便，dy 进行固定，0.0048 单位mm； ac固定pi/2。
% tx,ty,tz,rx,ry,rz,分别表示三个方向上的旋转和平移量。 处于inexI的5-10位置

ac=pi/2;
dy=0.0048;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % %操作仅仅为便于公式编写及读懂 % % % % % 

f=inexI(1);
dx=inexI(2);

u0=inexI(3);
v0=inexI(4);

% % % % %操作仅仅为便于公式编写及读懂 end % % % %

homouv=[uv 1]';

% Muvtoim=[dx 0 -dx*u0
%          0 dy -dy*v0
%          0 0 1];

% % % % % % % % % % % %完成表达式，包含ac % % % % % % % % % % % % 

Muvtoim=[dx dy*cos(ac) -dx*u0-dy*cos(ac)*v0
         0 dy*sin(ac) -dy*sin(ac)*v0
         0 0 1];     
     
homoxyi=Muvtoim*homouv;

xyi=[homoxyi(1) homoxyi(2)];


end


