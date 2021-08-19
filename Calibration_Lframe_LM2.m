function xf = Calibration_Lframe_LM2(x, xyzC, Lxyz)

% L型框架标定中使用的LM算法
% Calibration_LM2        去掉了没用的返回值、中间过程的输出以及DEBUG的内容
% Calibration_Lframe_LM2 L型框架标定中使用的LM算法，优化对象为L型框架上四个点在相机坐标系下的坐标 

% input:
% x     1*12     L型框架上四个点在相机坐标系下的坐标，按顺序排成一行. [ x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 ] 
% xyzC  4*3      L型框架上四个点在归一化平面上的坐标，一行表示一个点的坐标（z=1）. [ x y 1 ] 
% Lxyz  4*3      L型框架上四个点在世界坐标系下的坐标，一行表示一个点的坐标. [ x y z ] 

% ============ debug ================
% load input\inexIk0
% load input\uv_useful
% ===================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Default Options，初始化默认的参数

   MaxIter  = 0;       %   maximum number of iterations allowed，最大迭代次数
   ScaleD   = 1 ;        %   automatic scaling by D = diag(diag(J'*J))，单位矩阵I前面的系数
   FunTol   = 1e-7;      %   tolerace for final function value，目标函数的变化差足够小的值
   Lambda   = 0;         %   start with Newton iteration，Lambda的初始值
   Basdx    = 25e-9;     %   basic step dx for gradient evaluation，X的变化量足够小的值

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               INITIATION OF SOLUTION，再次初始化的过程
%               **********************
tic;
x = x(:) ;  %   Xo is a column vector
n = length(x);

bdx = Basdx*ones(n,1);      %   basic step dx for gradient evaluation
epsf  = FunTol(:);
maxit = MaxIter;    %   maximum permitted number of iterations
if maxit==0, maxit=100*n; end

l  = Lambda;
lc = 1;
r = Calibration_Lframe_LM2_lm4pt(x, xyzC, Lxyz); %   initial "residuals"，用x的初值计算出来的目标函数的分量，
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SS = r'*r;%目标函数

dx = zeros(n,1);
res= 1;
cnt=0;

[A,v] = getAv(x,r,bdx,xyzC, Lxyz) ; %调用的函数，计算出雅克比矩阵，然后用来计算A和v
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~

D = ScaleD(:);      %   CONSTANT SCALE CONTROL D
if isempty(D)
    D = diag(A);            %   automatic scaling
else
    ld = length(D);
    if ld==1
        D = abs(D)*ones(n,1); %   scalar of unique scaling
    elseif ld~=n
        error(['wrong number of scales D, lD = ',num2str(ld)])
    end
end
D(D<=0)=1;
T = sqrt(D);

Rlo=0.25;%R的范围
Rhi=0.75;

%               SOLUTION
%               ********    MAIN ITERATION CYCLE
while 1 %                   ********************

    cnt = cnt+1;%迭代次数
    
    d = diag(A);%A的对角线元素
    s = zeros(n,1); %长度与x相同的临时变量
    
%                           INTERNAL CYCLE
    while 1 %               ~~~~~~~~~~~~~~
        while 1
            UA = triu(A,1);%提取A的上三角矩阵
            A = UA'+UA+diag(d+l*D);%计算A+lambdaI，将A拆分为上三角和下三角，以及对角线，将I用D替换
            [U,p] = chol(A);        %   Choleski decomposition将A分解为U’*U，正定对称的即可，U是一个上三角矩阵
            %~~~~~~~~~~~~~~~
            if p==0, break, end
            l = 2*l;%l就是lambda
            if l==0, l=1; end
        end
        dx = U\(U'\v);              %   vector of x increments 公式中的delta
        vw = dx'*v;%异常情况处理，
        fin = -1;
        if vw<=0, break, end        %   The END

        for i=1:n    % 这个for循环主要实现的是A*dx的运算，z存储的是A*dx的第i个分量的值，这个A*dx的A是分解变化之前的、getAv计算出来的A.
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx; % dq = 2*v'*dx - dx'*A*dx ;
        s  = x-dx;
        rd = Calibration_Lframe_LM2_lm4pt(s, xyzC, Lxyz);
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        res = res+1;%重复寻找lambda的次数
        SSP = rd'*rd;%S'
        dS  = SS-SSP;
        fin = 1;
        if all((abs(dx)-bdx)<=0) || res>=maxit || abs(dS)<=epsf%限定迭代次数
            break                   %   The END
        end
        fin=0;
        if dS>=Rlo*dq, break, end
        A = U;
        y = .5;
        z = 2*vw-dS;
        if z>0, y=vw/z; end
        if y>.5, y=.5; end
        if y<.1, y=.1; end
        %计算lambdaC
        if l==0         % 计算 lambda 的边界值 lambda_c = 1/norm(inv(A))，作为 lambda 的新值，
                        % 这个inv(A)中的A是Choleski分解前的A.
                        % 一般来说，计算lambda_c的次数，在求解一个问题中不会超过2次.
                        % 由于原矩阵（Choleski分解前）A的对称性，整个计算过程中都将A表示为只有上三角的形式，
                        % 对A的操作也转化为对A的上三角部分做相应的操作，以节省计算次数，提高效率.
            y = 2*y;
            
            % 此时变量A的值为Choleski分解后上三角矩阵.后面两个for循环的作用是对上三角矩阵A求逆
            for i = 1:n % 先求对角线元素的逆
                A(i,i) = 1/A(i,i);
            end 
            for i = 2:n % 再求非对角线元素的逆
                ii = i-1;
                for j= 1:ii
                    A(j,i) = -A(j,j:ii)*A(j:ii,i).*A(i,i);
                end
            end
            
            % 得到A*A'（即Choleski分解前的A的逆）的上三角部分
            for i = 1:n 
                for j= i:n
                    A(i,j) = abs(A(i,j:n)*A(j,j:n)');
                end
            end
            
            % 计算Choleski分解前的A的逆的无穷范，即norm(inv(A),inf)
            l  = 0;
            tr = diag(A)'*D;
            for i = 1:n
                z = A(1:i,i)'*T(1:i)+z;
                if i<n
                    ii = i+1;
                    z  = A(i,ii:n)*T(ii:n)+z;
                end
                z = z*T(i);
                if z>l, l=z; end
            end
            if tr<l, l=tr; end
            l  = 1/l;
            lc = l;
        end
      %扩大lambda的值  
        l = l/y;
        if dS>0, dS=-1e300; break, end % dS>0 说明此次迭代的目标函数较上一轮有所减小，可以进入下一轮迭代. 
%                                        这里取dS=-1e300是为了不进入dS>Rhi*dq，因为这一步已经在dS<Rlo*dq的条件中了.
    end %  while            INTERNAL CYCLE LOOP
%                           ~~~~~~~~~~~~~~~~~~~

    if fin, break, end
    if dS>Rhi*dq
        l=l/2;
%         l = l / 10 ; %实际应用中，经常遇到l变到e+10量级，再慢慢除以2变到e+01量级的情况，加大分母的值可以加快下降的速度……
        if l<lc, l=0; end
    end
    SS=SSP;  x=s;  r=rd;
    [A,v] = getAv(x,r,bdx,xyzC, Lxyz) ;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % while                 MAIN ITERATION CYCLE LOOP
%                           *************************

if fin>0
    if dS>0
        SS = SSP;
        x  = s;
    end
end

xf = x;
if res>=maxit, cnt=-maxit; end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,v] = getAv(x,r,bdx,xyzC, Lxyz)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~  Calculate A, v, r

rc = r(:);
lx = length(x);
J  = zeros(length(r),lx);

parfor k = 1:lx
    dx = bdx(k);
    xd = x;
    xd(k) = xd(k)+dx;
    rd = Calibration_Lframe_LM2_lm4pt(xd, xyzC, Lxyz);
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    J(:,k)=((rd(:)-rc)/dx);

end

A = J'*J;
v = J'*r;

end % getAv
% --------------------------------------------------------------------



% end % LMFnlsq2