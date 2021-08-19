function Calibration_wand()

% ==== wand优化 ====
% 读入：内、外、畸变参数初值，wand在各镜头中的像素坐标
% 返回：优化后的内、外、畸变参数、M矩阵

camN = 16 ; % 总镜头数

% ==== Step1 读取wand在各镜头中的像素坐标 ====

% load input\wand %直接读取转换好的像素坐标
% frameN = zeros(1,camN) ;  %每个镜头的帧数
% for icam = 1:length(wand)
%     frameN(icam) = length(wand{icam}) ; 
% end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filename = 'CalWand5' ;
% datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
datapath ='input\20161018 展会数据\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
frameN = zeros(1,camN) ; %每个镜头的帧数

wand = [] ; %wand坐标，按镜头存储
fprintf('\n读取并转换wand数据：\n')
for icam = 1:camN
    wand{icam} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
    frameN(icam) = length(wand{icam}) ; 
    fprintf('第%d号镜头数据读取转换完毕\n',icam)
end

savepath = 'input\wand' ;
save(savepath,'wand')
fprintf('.vc转换后的wand文件已保存在input\\wand \n')


% ==== Step2 逐帧进行屏蔽杂点、wand筛选 ====

usecam = 1:camN ; %[1 4 5 8 11 12 15 16] ; % 正在使用的镜头编号
load input\shield ;              % shield,屏蔽区域
load input\inexIk0 ;             % inexIk,内外畸变参数初值，1*12camN

inexI = zeros(length(inexIk)/12,10) ;   %内外参，camN*10
for i = 1:length(inexIk)/12
    inexI(i,:) = inexIk(i*12-11:i*12-2) ;
end

usewand = wand(usecam) ;         % 正在使用的镜头数据
useinexI = inexI(usecam,:) ;     % 正在使用的镜头内外参
useshield = shield(usecam) ;     % 正在使用的镜头屏蔽信息

[ uv_useful, FrameUsed ] = Calibration_wand_Collect(usewand, useinexI, useshield) ; % 筛选


% ==== Step3 整体优化 ====

%parpool('local', 3);  %使用3核并行计算，4核会蓝屏……
%delete(gcp('nocreate'))
disp('LM优化开始');

% xf = Calibration_LM2(inexIk, uv_useful) ; %不输出中间过程
[xf, ~, ~, ~, ~] = Calibration_LM(inexIk, uv_useful) ; %输出中间过程

disp('LM优化结束');
save input\xf xf
disp('优化后的inexIk文件已保存在 input\xf ')

% load input\xf


% ==== Step4 调整优化后的坐标系位置到给定的L框架上 ====

% 拆分出畸变和内外参数
kk = [ xf(12*(1:camN)-1), xf(12*(1:camN)) ] ; % kk表示作用在实际坐标系上的畸变参数，与Cotex_setup中的意义不同，用于NK重建.
xf([12*(1:camN)-1; 12*(1:camN)]) = [] ; % 优化后的内外参

load input\Luv
load input\Lframexyz

inexI_hat = movedL(xf,Luv,xyz) ;  %参数调整


% ==== Step5 保存调整后的畸变参数、内参、外参、M矩阵 ====

% 保存畸变参数
dlmwrite('input\kk.txt',kk, 'delimiter','\t', 'precision', 10);
disp('畸变参数kk已保存在 input\kk.txt ')

% 保存内外参数
dlmwrite('input\inexI.txt',inexI_hat, 'delimiter','\t', 'precision', 10);
disp('内外参数inexI已保存在 input\inexI.txt ')

% 生成并保存NK用的M矩阵
M = [];
for i = 1:camN
    temp = buildM(inexI_hat(i,:)) ;
    M = [M temp'];
end
dlmwrite('input\M1.txt',M', 'delimiter','\t', 'precision', 10);
disp('M文件已保存在 input\M1.txt ')

% % 将畸变系数kk转换为Cotex的setup文件中的形式
% uv_useful_cam = cell(camN,1) ; % uv_useful是按帧存储的，将其转换为按镜头存储
% for iframe = 1:length(uv_useful)
%     uv = uv_useful{iframe} ;
%     for i = 1:size(uv,1)
%         uv_useful_cam{uv(i,1)} = [uv_useful_cam{uv(i,1)}, uv(i,2:5)] ;
%     end
% end
% k_setup = zeros(camN,2) ;
% for icam = 1:camN
%     uvr = uv_useful_cam{icam} ;
%     uvi = adddistortion(uvr, xf(icam*10-9:icam*10),  kk(icam,:)) ; %加畸变
%     k_setup(icam,:) = CptDistortion_uv(uvr,uvi,xf(icam*10-9:icam*10)) ;
% end
% dlmwrite('input\k_setup.txt',k_setup, 'delimiter','\t', 'precision', 10);
% disp('用于Cotex的setup文件中的畸变系数k已保存在 input\k_setup.txt ')
% 
% %生成Cotex的setup文件中的M矩阵
% M_setup = zeros(camN,11) ;
% for i = 1:size(inexI_hat,1)
%     M_setup(i,:) = NK2MAM(inexI_hat(i,:)) ;
% end
% dlmwrite('input\M_setup.txt',M_setup, 'delimiter','\t', 'precision', 10);
% disp('用于Cotex的setup文件中的M矩阵已保存在 input\M_setup.txt ')







