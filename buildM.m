function M = buildM(inexI)


% �������ܣ�����inexI���������������ֵ�����M����һ��1��10�С�
% % % % % ����M���󣬽��й�һ������% % % % % 

% % % �������������ʵ����ֵ����������ξ����M����
% inexI,��ʾ�ڲ��������1��4�У��� ��ʾ����������1��6�С�һ��һ��10�С�
% �ڲ�f dx u0 v0  ����inexI��1-4λ�á� Ϊʹ�÷��㣬dy ���й̶���0.0048 ��λmm�� ac�̶�pi/2��
% tx,ty,tz,rx,ry,rz,�ֱ��ʾ���������ϵ���ת��ƽ������ ����inexI��5-10λ��

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% fx,fy,ac,u0 v0,��ʾ�ڲε����������fx��fy��ʾ���࣬���ر�ʾ��ac��ֱ�ȣ����ȱ�ʾ
% in ��ʾ�ڲξ���3*3
% ex ��ʾ��ξ���3*4 [exR exT]
% exR �����ת��������
% exT ���ƽ�Ʋ�������
% M ��ʾM����,���Ÿ�ʽ1*11
% FM ��ʾ����M���� 3*4
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % %�˴���������ǻ��ȣ���Ҫ����һ��ת��

% % % % %��������Ϊ���ڹ�ʽ��д������ % % % % % 

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
rx=inexI(8)*pi/180; %�˴���������ǻ��ȣ�����һ��ת��
ry=inexI(9)*pi/180 ;%�˴���������ǻ��ȣ�����һ��ת��
rz=inexI(10)*pi/180 ;%�˴���������ǻ��ȣ�����һ��ת��
% % % % %��������Ϊ���ڹ�ʽ��д������ end % % % %

% % % % %�����ڲξ��󡣹�ʽ��˵�� % % % % %
% % �˴��ڲξ��󣬿��ܻ�����ֱ�����ò����������п�����Ҫ�Ľ�
in=[fx -fx/tan(ac) u0
    0 fy/sin(ac) v0  
    0 0 1]; 

% % % % %������ξ��󡣹�ʽ��˵�� % % % % % 



% % % XYZŷ��˳�� Rz * Rz *Rx 
% exR=[cos(rz)*cos(ry) sin(rz)*cos(rx)+cos(rz)*sin(ry)*sin(rx) sin(rz)*sin(rx)-cos(rz)*sin(ry)*cos(rx) 
%     -sin(rz)*cos(ry) cos(rz)*cos(rx)-sin(rz)*sin(ry)*sin(rx) cos(rz)*sin(rx)+sin(rz)*sin(ry)*cos(rx)
%     sin(ry) -cos(ry)*sin(rx) cos(ry)*cos(rx)];
% %     ex��ʽ�Ѿ��Ƶ���֤���������������Ѿ��Ƶ�������֤.

% % % ZYX ŷ��˳��....Rx* Rz *Rz 
exR=[cos(ry)*cos(rz) cos(ry)*sin(rz) -sin(ry)
    sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz) sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz) sin(rx)*cos(ry)
    cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz) cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz) cos(rx)*cos(ry)
    ]; % ex��ʽ�Ѿ��Ƶ���֤���������������Ѿ��Ƶ�������֤.  
% 
%   




TemexT=[tx;ty;tz];
exT=(-1)*exR*TemexT;

ex=[exR exT];   
 
FM=in*ex;
ss=FM/FM(3,4);
M=[ss(1,:) ss(2,:) ss(3,1) ss(3,2) ss(3,3)];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % ����������Գ���% % % % % % % 
% 

% r1=exR(:,1)
% r2=exR(:,2)
% r3=exR(:,3)
% r=[r1'*r1 r2'*r2 r3'*r3
%     r1'*r2 r2'*r3 r3'*r1]
% % % % % % % ����������Գ���end% % % % % % % 


end





