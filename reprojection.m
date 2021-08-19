
function puv = reprojection(M,xyz)
 

% 函数功能：输入3D坐标，和一个镜头的M矩阵，生成重投影坐标。
%          可同时计算多个点的重投影坐标
% 输入： xyz ,  3D坐标 n*3； 
%       M 矩阵，归一化之后的 1*11 或者 11*1 都可以
% 输出：puv，重投影坐标（uv形式）。n*2


% 单点的计算
    homoxyz=[xyz 1]; %homoxyz, 表示xyz的齐次坐标.横起来的形式。
    ss=[M(1) M(2) M(3) M(4);M(5) M(6) M(7) M(8);M(9) M(10) M(11) 1]*homoxyz'; %ss 临时性应用变量
    homopuv=ss'/ss(3); % %homopuv,表示重投影之后的2d坐标的齐次坐标，横起来的形式
    puv=[homopuv(1) homopuv(2)];  %puv,表示重投影之后的2d坐标,横着的形式。

% n = size(xyz,1) ; %点的个数
% homoxyz=[xyz, ones(n,1)]; % n*4, 表示xyz的齐次坐标.横起来的形式。
% ss=[M(1) M(2) M(3) M(4);M(5) M(6) M(7) M(8);M(9) M(10) M(11) 1]*homoxyz'; %ss 临时性应用变量，3行*n列
% homopuv=ss./repmat(ss(3,:),3,1); % %homopuv,表示重投影之后的2d坐标的齐次坐标，3*n
% puv = homopuv(1:2,:)';  %puv,表示重投影之后的2d坐标,横着的形式。

end
