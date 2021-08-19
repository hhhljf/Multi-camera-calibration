function Calibration_Lframe()

% L��ܱ궨

% clear

% ======================================================================
% ==== Step0. ���뾵ͷ���� ====
% ======================================================================

camN = 16 ; %ָ����ͷ����


% ======================================================================
% ==== Step1. L�Ϳ�ܱ궨 ====
% ���룺L�Ϳ���������꣬L�Ϳ���ڸ���ͷ�е��������꣬�ڲΣ����������
% ���أ��ڡ��⡢���������ֵ
% ======================================================================

% ==== Step1.1 ��ȡL�Ϳ�ܵ��������� ====

xyz=[0	0	0 
200	0	0
600	0	0
0	400	0
]; % �涨�õ�L�Ϳ����������
savepath = 'input\Lframexyz' ;
save(savepath,'xyz')
fprintf('L�Ϳ��3D����xyz�ļ��ѱ�����%s \n',savepath)

% ==== Step1.2 ��ȡL�Ϳ���ڸ���ͷ�е��������� ====

filename = 'CalFrame1' ;
% datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
datapath ='input\20161018 չ������\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
frameN = zeros(1,camN) ; %ÿ����ͷ��֡��
Luv = [] ; 
fprintf('\n��ȡ��ת��Lframe���ݣ�\n')
for icam = 1:camN
    Luv{icam} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
    frameN(icam) = length(Luv{icam}) ;
    fprintf('��%d�ž�ͷ���ݶ�ȡת�����\n',icam)
end
% �ɼ�������ֻ��60֡�����������趨

% ==== Step1.3 �����ӵ� ====

% ���L����ɼ�ʱ���ӵ�����shield������
shield = cell(1,camN) ; % cell:n*4��ÿ����ͷ���������б�Խǵ�����
shield_threshold = 20 ; %���ο�ĶԽ��߳���pix
for icam = 1:camN
    is = 1 ; %index_shield
    for iframe = 1:length(Luv{icam})
        while Luv{icam}{iframe}(2)<200  % v����С��200�ĵ㶼���ӵ㣬����ʵ������޸ģ���
            if iframe == 1
                point = Luv{icam}{iframe}(1:2) ;
                shield{icam}(is,:) = [point+shield_threshold*ones(1,2), point-shield_threshold*ones(1,2)] ;
                is = is+1 ;
            end
            Luv{icam}{iframe}(1:2) = [] ; %����
        end
    end
end
shield{7} = [740,40,642,0] ;%20161018 չ������ �ֶ��������ε�
savepath = ['input\',filename,'_shield'] ;
save(savepath,'Luv')
fprintf('�����ӵ���L������������ѱ����� %s \n', savepath)
savepath = ['input\','shield'] ;
save(savepath,'shield')
fprintf('���������ѱ����� %s \n', savepath)

% ==== Step1.4 ��ȡһ�����ڱ궨��L�Ϳ�ܵ�uv���� ====

% ��L��ܵ�uv���������ѡһ֡��Ϊ�궨�õ����ꡭ��
% Ҳ����ȡ��λ��������ȡƽ��ֵ�Ļ���Ҫȷ��û���ӵ㡭��
% iframe = randperm(length(Luv{1}),1) ;  %���ѡһ֡�������Ǻ��Ͻ�����
iframe = randperm(length(Luv{1})) ; 
iframe = iframe(1) ;%ʵ��60�������ȡһ֡
Luvnk = zeros(camN,8) ; % camN*8 ����ͷ��L���꣬û�����򣬷���Cotex_setup�еĸ�ʽ
for icam = 1:camN
    Luvnk(icam,:) = Luv{icam}{iframe} ;
end


% ==== Step1.5 �����ڲ� ====

%���ֲ��ڲΡ���
inI = zeros(camN,4) ;
inI(:,1) = 8      *ones(size(inI,1),1) ;
inI(:,2) = 0.0048 *ones(size(inI,1),1) ;
inI(:,3) = 960    *ones(size(inI,1),1) ;
inI(:,4) = 540    *ones(size(inI,1),1) ;


% ==== Step1.6 ����ͷ���μ������ ====

% ����ͷ��L����ת��Ϊ�궨�����ӿڵĸ�ʽ
Luv = zeros(4,2*camN) ; %����һ��u��v����Ĵ洢��ʽ
for i = 1:camN
    for k = 1:4
        Luv(k,2*(i-1)+1:2*i) = Luvnk(i,2*(k-1)+1:2*k) ;
    end
end
savepath = 'input\Luv' ;
save(savepath,'Luv')
fprintf('L�Ϳ��Luv�ļ��ѱ�����%s \n',savepath)

fprintf('\n��ʼ�����ⲿ������\n')
inexI = zeros(1,10*camN); % 1*10n ���������
for i = 1:camN
    inexI(i*10-9:i*10) = [ inI(i,:), LFrameCalibration_03(inI(i,:),Luv(:,i*2-1:i*2))] ; 
    fprintf('��%d�ž�ͷ��������������\n',i) ;
end


% ==== Step1.7 ��������ֵ����ֱ���ⲿ��������ֵ�� ====

% ������ƽ���ϵ�Luv��������
for i = 1:camN
    Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
end

% ������μ������ϵ��
c=2 ; %��������ĸ�����2��ʾ�����������Ϊk1,k2����
kk = zeros(camN,c) ; %����ϵ����һ�б�ʾһ����ͷ�Ļ���ϵ��
inexIk = zeros(1,camN*12) ; %����μӻ���ϵ��
for i = 1:camN
%     kk(i,:) = CptDistortion(Luv(:,i*2-1:i*2),xyz,inexI(i*10-9:i*10),c) ; %��������ֵ��Ч�����ã�����ֱ�Ӹ�ֵ
    kk(i,:) = [1.690005e-003, -3e-5] ; % ֱ���ⲿ��������ֵ
    inexIk(i*12-11:i*12) = [inexI(i*10-9:i*10), kk(i,:)] ;
end

savepath = 'input\inexIk0' ;
save(savepath,'inexIk')
fprintf('������������ֵinexIk�ļ��ѱ�����%s \n',savepath)








