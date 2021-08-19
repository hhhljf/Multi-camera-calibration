function d = Calibration_LM_wand_J(inexIk, uv_useful, uvi_useful, M0,  indparameter)

% LM�Ż������Ĺ�����ƽ��Ŀ�꺯�������������
% ����Jacobian����ר�á���
% �� Calibration_LM_wand ��ͬ���ǣ�ֻ������������ı�ľ�ͷ��M������������꣬������Ч��.
% 
% ��֡���㣬�ؽ�3D���꣬����˳���
% �ؽ���3D������ͶӰ��������ͷ������u,v�����ϵ����.
% �����䣬���䷽��ʵ������ -(����)-> ��������
% u,v �����ϵ����������������㣬��������Ϊ��
% ʵ������ -(����)-> �������� -(�ؽ�)-> 3D���� -(ͶӰ)-> ��������'
%                                  �Ƚ�3D�˳����      �ȽϹ�һ��ƽ���ϵ�2D���

% input:
% inexIk        1*12camN    [inexI, k]������μӻ���ϵ�� [ f dx u0 v0 tx ty tz rx ry rz k1 k2 ] 
% uv_useful     frameN*1    cell ÿ��cell�д洢һ֡��ͨ���ľ�ͷ�ż�a,b�����uv���꣬
%                           ÿ�и�ʽ[ ��ͷ��, ua, va, ub, vb ]
% uvi_useful    frameN*1    �����ı�ǰ������У�����������������꣬��uv_useful��Ӧ
% indparameter    ����      index parameter ����Jacobian����ʱ�����ı�Ĳ������±�
%                           -1 ��ʾû�в��������ı�
% M0            camN*11      �����ı�ǰ������ͷ��M����һ�б�ʾһ����ͷ�Ĳ���
% 
% output:
% d             1��          2D��3D�������

parameterN = 12 ; %ÿ����ͷ�Ĳ�������

if size(inexIk,2)==1, inexIk=inexIk(:); end

% ���㷢���ı�����ĸ���ͷ���ĸ�����
dpar  = rem(indparameter,parameterN) ; %�����ı���ǵڼ�������
if dpar==0
    dpar  = parameterN ;
    dcam = fix(indparameter/parameterN) ; %�����ı�ľ�ͷ��
else
    dcam = fix(indparameter/parameterN)+1 ;
end


frameN = length(uv_useful) ; %��֡��
camN = length(inexIk)/12 ;%camN ��ʾ��ͷ�ĸ���

inexIn = zeros(1,camN*10) ; % camN����ͷ�������
kn = zeros(1,camN*2) ; %camN����ͷ�Ļ���ϵ��
for i = 1:camN
    inexIn(10*i-9:10*i) = inexIk(12*i-11:12*i-2) ;
    kn(2*i-1:2*i) = inexIk(12*i-1:12*i) ;
end

% �������ı��������Σ�����Ҫ���¼���þ�ͷ��M����
if dpar<11
    M0(dcam,:) = buildM(inexIn(10*dcam-9:10*dcam)) ; % ��inexIn����M
end

% time_dis = 0 ;          %���������ʱ
% time_rebuild = 0 ;      %�ؽ���ʱ
% time_reprojection = 0 ; %��ͶӰ��ʱ

% ��֡���㣬��ͶӰ������uv��
d = []; %zeros( frameN*(1+), 1 ) ;
% d = zeros(frameN,1) ; %���Կ��ú��ء���
for iframe = 1:frameN

% tic    
    % ����У��
    uvi = uvi_useful{iframe} ; % uv_ideal ȥ�������������������
    M = M0( uvi(:,1), :) ;
    
    if dpar>10 %����ǻ�����������ı䣬����Ҫ���¼���þ�ͷ����������
        temp = find(uvi(:,1)==dcam) ;
        if ~isempty(temp)
            uv = uv_useful{iframe} ; % �������ʵ����������
            uvi(temp,2:5) = adddistortion(uv(temp,2:5), inexIn(uv(temp,1)*10-9:uv(temp,1)*10),  kn(uv(temp,1)*2-1:uv(temp,1)*2)) ; %�ӻ���
        end
    end
% time_dis = time_dis + toc ;    

% tic
    % �����ؽ�ʵ����̫���ˡ���ֱ�������о�ͷ��С���˰ɡ���
    tempuv = uvi(:,2:3)' ;
    tempuv = tempuv(:)' ;
    mxyza = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    tempuv = uvi(:,4:5)' ;
    tempuv = tempuv(:)' ;
    mxyzb = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    d = [d; 500-norm(mxyza-mxyzb)] ;
%    d(iframe) = abs(500-norm(mxyza-mxyzb)) ; 
%     d(iframe*2-1) = 500-norm(mxyza-mxyzb) ;  %���Կ��ú��ء���
% time_rebuild = time_rebuild + toc ;   

    % ��ͶӰ
% tic
    for i = 1:size(uvi,1)
        uva = reprojection( M(i,:), mxyza ) ;
        uvb = reprojection( M(i,:), mxyzb ) ;
        d = [d; (uvi(i,2:5)-[uva,uvb])'   ] ;
%         d = [d; norm(uvi(i,2:3)-uva); norm(uvi(i,4:5)-uvb) ] ; % d49569
%         d(iframe) = d(iframe) + sum(abs((uvi(i,2:5)-[uva,uvb]))) ; % d2981
    end
%     time_reprojection = time_reprojection + toc ;
end

% d=sum(abs(d)); %���Կ��ú��ء���

% fprintf('����У����ʱ��%f \n',time_dis)
% fprintf('�ؽ�Ӱ��ʱ��%f \n',time_rebuild)
% fprintf('��ͶӰ��ʱ��%f \n',time_reprojection)


end
