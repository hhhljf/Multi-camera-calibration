function uvr = adddistortion(uvi, inexI, k)

% ����ͷ�����������������uvi����ӻ��䣬�õ�ʵ���ĵ�����������uvr
% input��
% uvi   1*2n    uv_ideal������У����n����������������꣬ÿ����Ϊһ�����[u,v]���꣬[u1 v1 u2 v2 ...]
% inexI 1*10    ������ͷ�������������ʽΪ[ f dx u0 v0 tx ty tz rx ry rz ]
% k     1*c     c���������
% output:
% uvr   1*2n    uv_real��ʵ�ʲɼ�����n������������꣬ÿ����Ϊһ�����[u,v]���꣬[u1 v1 u2 v2 ...]

np = length(uvi)/2 ; %n_point����ĸ���
if length(k) ~= 2; error; end % ����ϵ���ĸ����ݶ�Ϊ����.
u0 = inexI(3) ;
v0 = inexI(4) ;

uvr = zeros(1,2*np) ; % ���ϻ�������������
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

% �������ܣ�����uv���ڲ���������ֵ�����ͼ����������xyi����
% uv,xyi 1*2

% inI,��ʾ�ڲ��������1��4�У�
% �ڲ�f dx u0 v0  ����inexI��1-4λ�á� Ϊʹ�÷��㣬dy ���й̶���0.0048 ��λmm�� ac�̶�pi/2��
% tx,ty,tz,rx,ry,rz,�ֱ��ʾ���������ϵ���ת��ƽ������ ����inexI��5-10λ��

ac=pi/2;
dy=0.0048;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % %��������Ϊ���ڹ�ʽ��д������ % % % % % 

f=inexI(1);
dx=inexI(2);

u0=inexI(3);
v0=inexI(4);

% % % % %��������Ϊ���ڹ�ʽ��д������ end % % % %

homouv=[uv 1]';

% Muvtoim=[dx 0 -dx*u0
%          0 dy -dy*v0
%          0 0 1];

% % % % % % % % % % % %��ɱ��ʽ������ac % % % % % % % % % % % % 

Muvtoim=[dx dy*cos(ac) -dx*u0-dy*cos(ac)*v0
         0 dy*sin(ac) -dy*sin(ac)*v0
         0 0 1];     
     
homoxyi=Muvtoim*homouv;

xyi=[homoxyi(1) homoxyi(2)];


end


