function Calibration_wand()

% ==== wand�Ż� ====
% ���룺�ڡ��⡢���������ֵ��wand�ڸ���ͷ�е���������
% ���أ��Ż�����ڡ��⡢���������M����

camN = 16 ; % �ܾ�ͷ��

% ==== Step1 ��ȡwand�ڸ���ͷ�е��������� ====

% load input\wand %ֱ�Ӷ�ȡת���õ���������
% frameN = zeros(1,camN) ;  %ÿ����ͷ��֡��
% for icam = 1:length(wand)
%     frameN(icam) = length(wand{icam}) ; 
% end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename = 'CalWand5' ;
% datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
datapath ='input\20161018 չ������\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
frameN = zeros(1,camN) ; %ÿ����ͷ��֡��

wand = [] ; %wand���꣬����ͷ�洢
fprintf('\n��ȡ��ת��wand���ݣ�\n')
for icam = 1:camN
    wand{icam} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
    frameN(icam) = length(wand{icam}) ; 
    fprintf('��%d�ž�ͷ���ݶ�ȡת�����\n',icam)
end

savepath = 'input\wand' ;
save(savepath,'wand')
fprintf('.vcת�����wand�ļ��ѱ�����input\\wand \n')


% ==== Step2 ��֡���������ӵ㡢wandɸѡ ====

usecam = 1:camN ; %[1 4 5 8 11 12 15 16] ; % ����ʹ�õľ�ͷ���
load input\shield ;              % shield,��������
load input\inexIk0 ;             % inexIk,������������ֵ��1*12camN

inexI = zeros(length(inexIk)/12,10) ;   %����Σ�camN*10
for i = 1:length(inexIk)/12
    inexI(i,:) = inexIk(i*12-11:i*12-2) ;
end

usewand = wand(usecam) ;         % ����ʹ�õľ�ͷ����
useinexI = inexI(usecam,:) ;     % ����ʹ�õľ�ͷ�����
useshield = shield(usecam) ;     % ����ʹ�õľ�ͷ������Ϣ

[ uv_useful, FrameUsed ] = Calibration_wand_Collect(usewand, useinexI, useshield) ; % ɸѡ


% ==== Step3 �����Ż� ====

%parpool('local', 3);  %ʹ��3�˲��м��㣬4�˻���������
%delete(gcp('nocreate'))
disp('LM�Ż���ʼ');

% xf = Calibration_LM2(inexIk, uv_useful) ; %������м����
[xf, ~, ~, ~, ~] = Calibration_LM(inexIk, uv_useful) ; %����м����

disp('LM�Ż�����');
save input\xf xf
disp('�Ż����inexIk�ļ��ѱ����� input\xf ')

% load input\xf


% ==== Step4 �����Ż��������ϵλ�õ�������L����� ====

% ��ֳ�������������
kk = [ xf(12*(1:camN)-1), xf(12*(1:camN)) ] ; % kk��ʾ������ʵ������ϵ�ϵĻ����������Cotex_setup�е����岻ͬ������NK�ؽ�.
xf([12*(1:camN)-1; 12*(1:camN)]) = [] ; % �Ż���������

load input\Luv
load input\Lframexyz

inexI_hat = movedL(xf,Luv,xyz) ;  %��������


% ==== Step5 ���������Ļ���������ڲΡ���Ρ�M���� ====

% ����������
dlmwrite('input\kk.txt',kk, 'delimiter','\t', 'precision', 10);
disp('�������kk�ѱ����� input\kk.txt ')

% �����������
dlmwrite('input\inexI.txt',inexI_hat, 'delimiter','\t', 'precision', 10);
disp('�������inexI�ѱ����� input\inexI.txt ')

% ���ɲ�����NK�õ�M����
M = [];
for i = 1:camN
    temp = buildM(inexI_hat(i,:)) ;
    M = [M temp'];
end
dlmwrite('input\M1.txt',M', 'delimiter','\t', 'precision', 10);
disp('M�ļ��ѱ����� input\M1.txt ')

% % ������ϵ��kkת��ΪCotex��setup�ļ��е���ʽ
% uv_useful_cam = cell(camN,1) ; % uv_useful�ǰ�֡�洢�ģ�����ת��Ϊ����ͷ�洢
% for iframe = 1:length(uv_useful)
%     uv = uv_useful{iframe} ;
%     for i = 1:size(uv,1)
%         uv_useful_cam{uv(i,1)} = [uv_useful_cam{uv(i,1)}, uv(i,2:5)] ;
%     end
% end
% k_setup = zeros(camN,2) ;
% for icam = 1:camN
%     uvr = uv_useful_cam{icam} ;
%     uvi = adddistortion(uvr, xf(icam*10-9:icam*10),  kk(icam,:)) ; %�ӻ���
%     k_setup(icam,:) = CptDistortion_uv(uvr,uvi,xf(icam*10-9:icam*10)) ;
% end
% dlmwrite('input\k_setup.txt',k_setup, 'delimiter','\t', 'precision', 10);
% disp('����Cotex��setup�ļ��еĻ���ϵ��k�ѱ����� input\k_setup.txt ')
% 
% %����Cotex��setup�ļ��е�M����
% M_setup = zeros(camN,11) ;
% for i = 1:size(inexI_hat,1)
%     M_setup(i,:) = NK2MAM(inexI_hat(i,:)) ;
% end
% dlmwrite('input\M_setup.txt',M_setup, 'delimiter','\t', 'precision', 10);
% disp('����Cotex��setup�ļ��е�M�����ѱ����� input\M_setup.txt ')







