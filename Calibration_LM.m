function [xf, SS, cnt, res, XY] = Calibration_LM(inexIk, uv_useful)

% 优化内外畸变参数的LM算法

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

   Display  = [1,0];     %   no print of iterations
   MaxIter  = 300;       %   maximum number of iterations allowed
   ScaleD   = abs(1./inexIk)' ;        %   automatic scaling by D = diag(diag(J'*J))
   FunTol   = 1e-7;      %   tolerace for final function value
   Printf   = 'printit'; %   disply intermediate results
   Trace    = 0;         %   don't save  intermediate results
   Lambda   = 1;         %   start with Newton iteration
   Basdx    = 25e-9;     %   basic step dx for gradient evaluation

% FUN = 'wandLM2_dis2' ;
%   [Xf,Ssq,CNT,Res,XY] = LMFnlsq2(FUN,Xo,Options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               INITIATION OF SOLUTION
%               **********************
tic;
x = inexIk(:) ;  %   Xo is a column vector
n = length(x);

bdx = abs(x*Basdx);        %   basic step dx for gradient evaluation
% lb  = length(bdx);
% if lb==1
%     bdx = bdx*ones(n,1);
% elseif lb~=n
%     error(['Dimensions of vector dx ',num2str(lb),'~=',num2str(n)]);
% end

epsf  = FunTol(:);
ipr   = Display;
maxit = MaxIter;    %   maximum permitted number of iterations
if maxit==0, maxit=100*n; end
printf= Printf;

l  = Lambda;
lc = 1;
[r, M0, uvi_useful ] = Calibration_LM_wand(x, uv_useful); %   initial "residuals"
% r  = feval(FUN,x);          %   initial "residuals"
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SS = r'*r;

feval(printf,ipr,-1);       %   Table header
dx = zeros(n,1);
res= 1;
cnt=0;
feval(printf,ipr,cnt,res,SS,x,dx,l,lc) %    Initial state

[A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0,ipr) ; 
% [A,v] = getAv(FUN,JAC,x,r,bdx,ipr);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trcXY = Trace;      %   iteration tracing
if trcXY
    XY = zeros(n,maxit);
    XY(:,1) = x;
else
    XY = [];
end

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
    if cnt>0
        feval(printf,ipr,cnt,res,SS,x,dx,l,lc)
    end
    cnt = cnt+1;
    if trcXY, XY(:,cnt+1)=x; end
    
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

        for i=1:n    % 这个for循环主要实现的是A*dx的运算，z存储的是A*dx的第i个分量的值，这个A*dx的A是分解变化之前的、getAv计算出来的A.
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx; % dq = 2*v'*dx - dx'*A*dx ;
        s  = x-dx;
        [rd, M0, uvi_useful ] = Calibration_LM_wand(s, uv_useful);
%         rd = feval(FUN,s);
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
    [A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0,ipr) ;
%     [A,v] = getAv(FUN,JAC,x,r,bdx,ipr);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % while                 MAIN ITERATION CYCLE LOOP
%                           *************************

if fin>0
    if dS>0
        SS = SSP;
        x  = s;
    end
end
if ipr(1)~=0
    disp(' ');
    feval(printf,sign(ipr),cnt,res,SS,x,dx,l,lc)
end
xf = x;
if trcXY, XY(:,cnt+2)=x; end
XY(:,cnt+3:end) = [];
if res>=maxit, cnt=-maxit; end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,v] = getAv(x,r,bdx,uv_useful,uvi_useful, M0,ipr)
% function [A,v] = getAv(FUN,JAC,x,r,bdx,ipr)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Calculate A, v, r

rc = r(:);
lx = length(x);
J  = zeros(length(r),lx);
tic
for k = 1:lx
    dx = bdx(k);
    xd = x;
    xd(k) = xd(k)+dx;
    rd = Calibration_LM_wand_J(xd, uv_useful, uvi_useful, M0,  k);
%     [rd, nouse, nouse2 ] = Calibration_LM_wand(xd, uv_useful) ;
%     rd = feval(FUN,xd);
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    J(:,k)=((rd(:)-rc)/dx);
    
% output columns of Jacobian matrix   
  f =' %7.4f %7.4f %7.4f %7.4f';
  if ipr(1)~=0 && ipr(2)>0
    fprintf(' Column #%3d\n',k);
    tim  = toc;
    mins = floor(tim/60);
    secs = tim-mins*60;
    fprintf(' Elapsed time  =%4d min%5.1f sec\n',mins,secs)
    if ipr(2)==2
      fprintf([ f f f f '\n'], J(:,k)');
    end
    fprintf('\n');
  end

end
toc
A = J'*J;
v = J'*r;

end % getAv
% --------------------------------------------------------------------


function printit(ipr,cnt,res,SS,x,dx,l,lc)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Printing of intermediate results
% For length(ipr)==1: 
%  ipr(1) <  0  do not print lambda columns
%         =  0  do not print at all
%         >  0  print every (ipr)th iteration
%  ipr(2) =  0  do not prit time neither J
% For length(ipr)>1 && ipr(1)~=0:
%               print also elapsed time between displays
%  cnt = -1  print out the header
%        >0  print out results
if ipr(1)~=0
   if cnt<0                 %   table header
      disp('')
      nch = 50+(ipr(1)>0)*25;
      disp(char('*'*ones(1,nch)))
      fprintf('  itr  nfJ   SUM(r^2)        x           dx');
      if ipr(1)>0
          fprintf('           l           lc');
      end
      fprintf('\n');
      disp(char('*'*ones(1,nch)))
      disp('')
   else                     %   iteration output
      if rem(cnt,ipr(1))==0
          if ipr(2)>0
              tim  = toc;
              mins = floor(tim/60);
              secs = tim-mins*60;
              fprintf('\n Elapsed time  =%4d min%5.1f sec\n',mins,secs)
          end
          f='%12.4e ';
          if ipr(1)>0
             fprintf(['%5d %5d ' f f f f f '\n'],...
                 cnt,res,SS, x(1),dx(1),l,lc);
          else
             fprintf(['%5d %5d ' f f f '\n'],...
                 cnt,res,SS, x(1),dx(1));
          end
%           for k=2:length(x)
          for k=2:12  %只显示一个镜头的变量
             fprintf([blanks(25) f f '\n'],x(k),dx(k));
          end
      end
   end
end
end % printit
% end % LMFnlsq2