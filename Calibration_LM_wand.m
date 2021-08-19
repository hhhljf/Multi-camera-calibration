function [d, M0, uvi_useful ] = Calibration_LM_wand(inexIk, uv_useful)

% LM�Ż������Ĺ�����ƽ��Ŀ�꺯�������������
% ��֡���㣬�ؽ�3D���꣬����˳���
% ���ؽ���3D������ͶӰ��������ͷ������u,v�����ϵ����.
% �����䣬���䷽��ʵ������ -(����)-> ��������
% u,v �����ϵ����������������㣬��������Ϊ��
% ʵ������ -(����)-> �������� -(�ؽ�)-> 3D���� -(ͶӰ)-> ��������'
%                                  �Ƚ�3D�˳����      �ȽϹ�һ��ƽ���ϵ�2D���

% input:
% inexIk        1*12camN    [inexI, k]������μӻ���ϵ�� [ f dx u0 v0 tx ty tz rx ry rz k1 k2 ] 
% uv_useful     frameN*1    cell ÿ��cell�д洢һ֡��ͨ���ľ�ͷ�ż�a,b�����uv���꣬
%                           ÿ�и�ʽ[ ��ͷ��, ua, va, ub, vb ]
% indparameter    ����      index parameter ����Jacobian����ʱ�����ı�Ĳ������±�
%                           -1 ��ʾû�в��������ı�
% 
% output:
% d             1��          2D��3D�������
% M0            camN*11      ���μ�����M����
% uvi_useful     frameN*1     ����У�������������

parameterN = 12 ; %ÿ����ͷ�Ĳ�������
uvi_useful = cell(size(uv_useful)) ; %����У�������������

% if size(inexIk,2)==1, inexIk=inexIk(:); end

frameN = length(uv_useful) ; %��֡��
camN = length(inexIk)/parameterN ;%camN ��ʾ��ͷ�ĸ���

inexIn = zeros(1,camN*10) ; % camN����ͷ�������
kn = zeros(1,camN*2) ; %camN����ͷ�Ļ���ϵ��
for i = 1:camN
    inexIn(10*i-9:10*i) = inexIk(12*i-11:12*i-2) ;
    kn(2*i-1:2*i) = inexIk(12*i-1:12*i) ;
end

% ��inexIn����M
M0 = zeros(camN,11) ;
for i = 1:camN
    M0(i,:) = buildM(inexIn(10*i-9:10*i)) ;
end

% time_dis = 0 ;          %���������ʱ
% time_rebuild = 0 ;      %�ؽ���ʱ
% time_reprojection = 0 ; %��ͶӰ��ʱ

% ��֡���㣬��ͶӰ������uv��
d = []; %zeros( frameN*(1+), 1 ) ;
% d = zeros(frameN,1) ; %���Կ��ú��ء���

for iframe = 1:frameN
    uv = uv_useful{iframe} ;
    M = M0( uv(:,1), :) ;
% tic    
    % ����У��
    uvi = zeros(size(uv)) ; % uv_ideal ȥ�������������������
    uvi(:,1) = uv(:,1) ;
    for iuv = 1:size(uv,1) ;
        uvi(iuv,2:5) = adddistortion(uv(iuv,2:5), inexIn(uv(iuv,1)*10-9:uv(iuv,1)*10),  kn(uv(iuv,1)*2-1:uv(iuv,1)*2)) ; %�ӻ���
    end
    uvi_useful{iframe} = uvi ;
% time_dis = time_dis + toc ;    

%     usefulNum = size(uv,1) ;
%     % �����ؽ�������a,b���3D����
%     k = 1 ;
%     xyza = zeros(usefulNum*(usefulNum-1)/2, 3) ; % a���Ӧ�� [ x, y, z ]
%     xyzb = zeros(usefulNum*(usefulNum-1)/2, 3) ; % b���Ӧ�� [ x, y, z ]
%     for i = 1:usefulNum-1
%         for j = i+1:usefulNum
%             xyza(k,:) = rebulid_3D_UnspecCam_LM( [uv(i,2:3), uv(j,2:3)], M([i,j],:) )  ;
%             xyzb(k,:) = rebulid_3D_UnspecCam_LM( [uv(i,4:5), uv(j,4:5)], M([i,j],:) )  ;
%             k = k+1 ;
%         end
%     end
%     mxyza = mean(xyza,1) ;
%     mxyzb = mean(xyzb,1) ;
%     d = [d; 38.7298*(500-norm(mxyza-mxyzb))] ;

% tic
    % �����ؽ�ʵ����̫���ˡ���ֱ�������о�ͷ��С���˰ɡ���
    tempuv = uvi(:,2:3)' ;
    tempuv = tempuv(:)' ;
    mxyza = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    tempuv = uvi(:,4:5)' ;
    tempuv = tempuv(:)' ;
    mxyzb = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
% %     d = [d; 38.7298*(500-norm(mxyza-mxyzb))] ;
    d = [d; 500-norm(mxyza-mxyzb)] ;
%     d(iframe) = abs(500-norm(mxyza-mxyzb)) ; %���Կ��ú��ء���
% time_rebuild = time_rebuild + toc ;    

    % ��ͶӰ
% tic
    for i = 1:size(uv,1)
        uva = reprojection( M(i,:), mxyza ) ;
        uvb = reprojection( M(i,:), mxyzb ) ;
        d = [d; (uvi(i,2:5)-[uva,uvb])'   ] ;
%         d = [d; norm(uvi(i,2:3)-uva); norm(uvi(i,4:5)-uvb) ] ; % d49569
%         d(iframe) = d(iframe) + sum(abs((uvi(i,2:5)-[uva,uvb]))) ; % d2981
    end
    
%     d=sum(abs(d)); %���Կ��ú��ء���
    
% time_reprojection = time_reprojection + toc ;
end

% fprintf('����У����ʱ��%f \n',time_dis)
% fprintf('�ؽ�Ӱ��ʱ��%f \n',time_rebuild)
% fprintf('��ͶӰ��ʱ��%f \n',time_reprojection)

% ------------------------------------------------
% % ����L��ܵ���Ϣ
% 
% xyz=[0	0	0 
% 200	0	0
% 600	0	0
% 0	400	0
% ]; % �涨�õ�L�Ϳ����������
% 
% Luv = load('C:\Users\Boat\Desktop\ʵ��\data_uvc\CalFrame1_shield') ; %ֱ�Ӷ�ȡת���õ�L����
% % Luv = load('input\Luv') ;
% Luv = Luv.uv ;
% iframe = randperm(length(Luv{1}),1) ;
% % iframe = randperm(min(frameN),1) ;  %���ѡһ֡������
% Luvr = zeros(camN,8) ; % Luv_real ����ͷʵ���ĵ���L����
% for icam = 1:camN
%     Luvr(icam,:) = Luv{icam}{iframe} ;
% end
% 
% % �����ؽ���3D���
% Luv = zeros(4,2*camN) ; %����ͷ��L����ת��Ϊ�궨�����ӿڵĸ�ʽ
% for i = 1:camN
%     for k = 1:4
%         Luv(k,2*(i-1)+1:2*i) = Luvr(i,2*(k-1)+1:2*k) ;
%     end
%     Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
% end
% 
% dd = zeros(4,1) ; % L�Ϳ���ĸ����3D���
% for k = 1:4
%     dd(k) = norm( xyz(k,:) - rebulid_3D_UnspecCam_LM(Luv(k,:),M0) ) ;
% end
% 
% d = [d; sqrt(iframe/15)*dd] ;

% ������ͶӰ��2D���
% Luvi = zeros(camN,8) ; % Luv_ideal ����ͷ��ͶӰ��L����
% for icam = 1:camN
%     for i = 1:4
%         Luvi(icam,i*2-1:i*2) = reprojection(M0(icam,:),xyz(i,:)) ;
%     end
% end
% 
% Luvr = Luvr(:) ;
% Luvi = Luvi(:) ;
% 
% d = [d; (Luvr - Luvi)] ;


