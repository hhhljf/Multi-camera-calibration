function xf = Calibration_Lframe_LM2(x, xyzC, Lxyz)

% L�Ϳ�ܱ궨��ʹ�õ�LM�㷨
% Calibration_LM2        ȥ����û�õķ���ֵ���м���̵�����Լ�DEBUG������
% Calibration_Lframe_LM2 L�Ϳ�ܱ궨��ʹ�õ�LM�㷨���Ż�����ΪL�Ϳ�����ĸ������������ϵ�µ����� 

% input:
% x     1*12     L�Ϳ�����ĸ������������ϵ�µ����꣬��˳���ų�һ��. [ x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 ] 
% xyzC  4*3      L�Ϳ�����ĸ����ڹ�һ��ƽ���ϵ����꣬һ�б�ʾһ��������꣨z=1��. [ x y 1 ] 
% Lxyz  4*3      L�Ϳ�����ĸ�������������ϵ�µ����꣬һ�б�ʾһ���������. [ x y z ] 

% ============ debug ================
% load input\inexIk0
% load input\uv_useful
% ===================================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Default Options����ʼ��Ĭ�ϵĲ���

   MaxIter  = 0;       %   maximum number of iterations allowed������������
   ScaleD   = 1 ;        %   automatic scaling by D = diag(diag(J'*J))����λ����Iǰ���ϵ��
   FunTol   = 1e-7;      %   tolerace for final function value��Ŀ�꺯���ı仯���㹻С��ֵ
   Lambda   = 0;         %   start with Newton iteration��Lambda�ĳ�ʼֵ
   Basdx    = 25e-9;     %   basic step dx for gradient evaluation��X�ı仯���㹻С��ֵ

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%               INITIATION OF SOLUTION���ٴγ�ʼ���Ĺ���
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
r = Calibration_Lframe_LM2_lm4pt(x, xyzC, Lxyz); %   initial "residuals"����x�ĳ�ֵ���������Ŀ�꺯���ķ�����
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SS = r'*r;%Ŀ�꺯��

dx = zeros(n,1);
res= 1;
cnt=0;

[A,v] = getAv(x,r,bdx,xyzC, Lxyz) ; %���õĺ�����������ſ˱Ⱦ���Ȼ����������A��v
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

Rlo=0.25;%R�ķ�Χ
Rhi=0.75;

%               SOLUTION
%               ********    MAIN ITERATION CYCLE
while 1 %                   ********************

    cnt = cnt+1;%��������
    
    d = diag(A);%A�ĶԽ���Ԫ��
    s = zeros(n,1); %������x��ͬ����ʱ����
    
%                           INTERNAL CYCLE
    while 1 %               ~~~~~~~~~~~~~~
        while 1
            UA = triu(A,1);%��ȡA�������Ǿ���
            A = UA'+UA+diag(d+l*D);%����A+lambdaI����A���Ϊ�����Ǻ������ǣ��Լ��Խ��ߣ���I��D�滻
            [U,p] = chol(A);        %   Choleski decomposition��A�ֽ�ΪU��*U�������ԳƵļ��ɣ�U��һ�������Ǿ���
            %~~~~~~~~~~~~~~~
            if p==0, break, end
            l = 2*l;%l����lambda
            if l==0, l=1; end
        end
        dx = U\(U'\v);              %   vector of x increments ��ʽ�е�delta
        vw = dx'*v;%�쳣�������
        fin = -1;
        if vw<=0, break, end        %   The END

        for i=1:n    % ���forѭ����Ҫʵ�ֵ���A*dx�����㣬z�洢����A*dx�ĵ�i��������ֵ�����A*dx��A�Ƿֽ�仯֮ǰ�ġ�getAv���������A.
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx; % dq = 2*v'*dx - dx'*A*dx ;
        s  = x-dx;
        rd = Calibration_Lframe_LM2_lm4pt(s, xyzC, Lxyz);
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        res = res+1;%�ظ�Ѱ��lambda�Ĵ���
        SSP = rd'*rd;%S'
        dS  = SS-SSP;
        fin = 1;
        if all((abs(dx)-bdx)<=0) || res>=maxit || abs(dS)<=epsf%�޶���������
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
        %����lambdaC
        if l==0         % ���� lambda �ı߽�ֵ lambda_c = 1/norm(inv(A))����Ϊ lambda ����ֵ��
                        % ���inv(A)�е�A��Choleski�ֽ�ǰ��A.
                        % һ����˵������lambda_c�Ĵ����������һ�������в��ᳬ��2��.
                        % ����ԭ����Choleski�ֽ�ǰ��A�ĶԳ��ԣ�������������ж���A��ʾΪֻ�������ǵ���ʽ��
                        % ��A�Ĳ���Ҳת��Ϊ��A�������ǲ�������Ӧ�Ĳ������Խ�ʡ������������Ч��.
            y = 2*y;
            
            % ��ʱ����A��ֵΪCholeski�ֽ�������Ǿ���.��������forѭ���������Ƕ������Ǿ���A����
            for i = 1:n % ����Խ���Ԫ�ص���
                A(i,i) = 1/A(i,i);
            end 
            for i = 2:n % ����ǶԽ���Ԫ�ص���
                ii = i-1;
                for j= 1:ii
                    A(j,i) = -A(j,j:ii)*A(j:ii,i).*A(i,i);
                end
            end
            
            % �õ�A*A'����Choleski�ֽ�ǰ��A���棩�������ǲ���
            for i = 1:n 
                for j= i:n
                    A(i,j) = abs(A(i,j:n)*A(j,j:n)');
                end
            end
            
            % ����Choleski�ֽ�ǰ��A������������norm(inv(A),inf)
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
      %����lambda��ֵ  
        l = l/y;
        if dS>0, dS=-1e300; break, end % dS>0 ˵���˴ε�����Ŀ�꺯������һ��������С�����Խ�����һ�ֵ���. 
%                                        ����ȡdS=-1e300��Ϊ�˲�����dS>Rhi*dq����Ϊ��һ���Ѿ���dS<Rlo*dq����������.
    end %  while            INTERNAL CYCLE LOOP
%                           ~~~~~~~~~~~~~~~~~~~

    if fin, break, end
    if dS>Rhi*dq
        l=l/2;
%         l = l / 10 ; %ʵ��Ӧ���У���������l�䵽e+10����������������2�䵽e+01������������Ӵ��ĸ��ֵ���Լӿ��½����ٶȡ���
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