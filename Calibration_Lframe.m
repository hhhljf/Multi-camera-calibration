function Calibration_Lframe()

% L框架标定

% clear

% ======================================================================
% ==== Step0. 读入镜头个数 ====
% ======================================================================

camN = 16 ; %指定镜头数量


% ======================================================================
% ==== Step1. L型框架标定 ====
% 读入：L型框架世界坐标，L型框架在各镜头中的像素坐标，内参，畸变参数；
% 返回：内、外、畸变参数初值
% ======================================================================

% ==== Step1.1 读取L型框架的世界坐标 ====

xyz=[0	0	0 
200	0	0
600	0	0
0	400	0
]; % 规定好的L型框架世界坐标
savepath = 'input\Lframexyz' ;
save(savepath,'xyz')
fprintf('L型框架3D坐标xyz文件已保存在%s \n',savepath)

% ==== Step1.2 读取L型框架在各镜头中的像素坐标 ====

filename = 'CalFrame1' ;
% datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
datapath ='input\20161018 展会数据\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
frameN = zeros(1,camN) ; %每个镜头的帧数
Luv = [] ; 
fprintf('\n读取并转换Lframe数据：\n')
for icam = 1:camN
    Luv{icam} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
    frameN(icam) = length(Luv{icam}) ;
    fprintf('第%d号镜头数据读取转换完毕\n',icam)
end
% 采集的数据只有60帧，可以自由设定

% ==== Step1.3 屏蔽杂点 ====

% 标记L坐标采集时的杂点区域shield并屏蔽
shield = cell(1,camN) ; % cell:n*4，每个镜头屏蔽区域的斜对角点坐标
shield_threshold = 20 ; %屏蔽框的对角线长度pix
for icam = 1:camN
    is = 1 ; %index_shield
    for iframe = 1:length(Luv{icam})
        while Luv{icam}{iframe}(2)<200  % v坐标小于200的点都是杂点，根据实际情况修改！！
            if iframe == 1
                point = Luv{icam}{iframe}(1:2) ;
                shield{icam}(is,:) = [point+shield_threshold*ones(1,2), point-shield_threshold*ones(1,2)] ;
                is = is+1 ;
            end
            Luv{icam}{iframe}(1:2) = [] ; %屏蔽
        end
    end
end
shield{7} = [740,40,642,0] ;%20161018 展会数据 手动补的屏蔽点
savepath = ['input\',filename,'_shield'] ;
save(savepath,'Luv')
fprintf('屏蔽杂点后的L框架像素坐标已保存在 %s \n', savepath)
savepath = ['input\','shield'] ;
save(savepath,'shield')
fprintf('屏蔽区域已保存在 %s \n', savepath)

% ==== Step1.4 提取一组用于标定的L型框架的uv坐标 ====

% 在L框架的uv坐标中随便选一帧作为标定用的坐标……
% 也可以取中位数……，取平均值的话需要确保没有杂点……
% iframe = randperm(length(Luv{1}),1) ;  %随便选一帧……不是很严谨……
iframe = randperm(length(Luv{1})) ; 
iframe = iframe(1) ;%实现60里面随便取一帧
Luvnk = zeros(camN,8) ; % camN*8 各镜头的L坐标，没有排序，仿照Cotex_setup中的格式
for icam = 1:camN
    Luvnk(icam,:) = Luv{icam}{iframe} ;
end


% ==== Step1.5 读入内参 ====

%用手册内参……
inI = zeros(camN,4) ;
inI(:,1) = 8      *ones(size(inI,1),1) ;
inI(:,2) = 0.0048 *ones(size(inI,1),1) ;
inI(:,3) = 960    *ones(size(inI,1),1) ;
inI(:,4) = 540    *ones(size(inI,1),1) ;


% ==== Step1.6 按镜头依次计算外参 ====

% 各镜头的L坐标转换为标定函数接口的格式
Luv = zeros(4,2*camN) ; %换了一种u和v坐标的存储方式
for i = 1:camN
    for k = 1:4
        Luv(k,2*(i-1)+1:2*i) = Luvnk(i,2*(k-1)+1:2*k) ;
    end
end
savepath = 'input\Luv' ;
save(savepath,'Luv')
fprintf('L型框架Luv文件已保存在%s \n',savepath)

fprintf('\n开始计算外部参数：\n')
inexI = zeros(1,10*camN); % 1*10n 内外参向量
for i = 1:camN
    inexI(i*10-9:i*10) = [ inI(i,:), LFrameCalibration_03(inI(i,:),Luv(:,i*2-1:i*2))] ; 
    fprintf('第%d号镜头内外参数计算完毕\n',i) ;
end


% ==== Step1.7 计算畸变初值（或直接外部读入畸变初值） ====

% 对像素平面上的Luv坐标排序
for i = 1:camN
    Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
end

% 用内外参计算畸变系数
c=2 ; %畸变参数的个数，2表示畸变参数个数为k1,k2两个
kk = zeros(camN,c) ; %畸变系数，一行表示一个镜头的畸变系数
inexIk = zeros(1,camN*12) ; %内外参加畸变系数
for i = 1:camN
%     kk(i,:) = CptDistortion(Luv(:,i*2-1:i*2),xyz,inexI(i*10-9:i*10),c) ; %计算畸变初值，效果不好，建议直接赋值
    kk(i,:) = [1.690005e-003, -3e-5] ; % 直接外部读入畸变初值
    inexIk(i*12-11:i*12) = [inexI(i*10-9:i*10), kk(i,:)] ;
end

savepath = 'input\inexIk0' ;
save(savepath,'inexIk')
fprintf('内外畸变参数初值inexIk文件已保存在%s \n',savepath)








