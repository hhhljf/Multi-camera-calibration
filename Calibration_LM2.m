function xf = Calibration_LM2(inexIk, uv_useful)

% 优化内外畸变参数的LM算法
% Calibration_LM2 去掉了没用的返回值、中间过程的输出以及DEBUG的内容

% input:
% inexIk    1*12camN    各镜头的内外参数，按顺序排成一行，每12个为一个镜头的参数. [ f dx u0 v0 tx ty tz rx ry rz k1 k2 ] 
% uv_useful frameN*1 cell 每个cell中存储一帧中通过的镜头号及a,b两点的uv坐标，
%                         每行格式[ 镜头号, ua, va, ub, vb ]

% ============ debug ================
% load input\inexIk0
% load input\uv_useful
% ===================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Default Options

   MaxIter  = 300;       %   maximum number of iterations allowed
   ScaleD   = abs(1./inexIk)' ;        %   automatic scaling by D = diag(diag(J'*J))
   FunTol   = 1e-7;      %   tolerace for final function value
   Lambda   = 2;         %   start with Newton iteration
   Basdx    = 25e-5;     %   basic step dx for gradient evaluation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               INITIATION OF SOLUTION
%               **********************
tic;
x = inexIk(:) ;  %   Xo is a column vector
n = length(x);

bdx = abs(x*Basdx);        %   basic step dx for gradient evaluation
epsf  = FunTol(:);
maxit = MaxIter;    %   maximum permitted number of iterations
if maxit==0, maxit=100*n; end

l  = Lambda;
lc = 1;
[r, M0, uvi_useful ] = Calibration_LM_wand(x, uv_useful); %   initial "residuals"
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SS = r'*r;

dx = zeros(n,1);
res= 1;
cnt=0;

[A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0) ; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Rlo=0.25;
Rhi=0.75;

%               SOLUTION
%               ********    MAIN ITERATION CYCLE
while 1 %                   ********************

    cnt = cnt+1;
    
    d = diag(A);
    s = zeros(n,1); %长度与x相同的临时变量
    
%                           INTERNAL CYCLE
    while 1 %               ~~~~~~~~~~~~~~
        while 1
            UA = triu(A,1);
            A = UA'+UA+diag(d+l*D);
            [U,p] = chol(A);        %   Choleski decomposition
            %~~~~~~~~~~~~~~~
            if p==0, break, end
            l = 2*l;
            if l==0, l=1; end
        end
        dx = U\(U'\v);              %   vector of x increments
        vw = dx'*v;
        fin = -1;
        if vw<=0, break, end        %   The END

        for i=1:n    % 这个for循环主要实现的是A*dx的运算，z存储的是A*dx的第i个分量的值. 这个A*dx的A是分解变化之前的、getAv计算出来的A.
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx; % dq = 2*v'*dx - dx'*A*dx ;
        s  = x-dx;
        [rd, M0, uvi_useful ] = Calibration_LM_wand(s, uv_useful);
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        res = res+1;
        SSP = rd'*rd;
        dS  = SS-SSP;
        fin = 1;
        if all((abs(dx)-bdx)<=0) || res>=maxit || abs(dS)<=epsf
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
            
            % 计算A*A'（即Choleski分解前的A的逆）的上三角部分，含对角线
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
    [A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0) ;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

function [A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Calculate A, v, r

rc = r(:);
lx = length(x);
J  = zeros(length(r),lx);

parfor k = 1:lx
    dx = bdx(k);
    xd = x;
    xd(k) = xd(k)+dx;
    rd = Calibration_LM_wand_J(xd, uv_useful, uvi_useful, M0,  k);
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    J(:,k)=((rd(:)-rc)/dx);

end

A = J'*J;
v = J'*r;

end % getAv
% --------------------------------------------------------------------

% end % LMFnlsq2