function [ uv_useful, FrameUsed ] = Calibration_wand_Collect(wand, inexI, shield)

% wand标定中的筛选部分
% 主要包含屏蔽杂点、第一步筛选和第二步筛选，共三个步骤
% 
% input:
% wand      1*camN   cell型，各镜头的wand像素坐标，三层cell结构……
%                    每个cell: 1*frameN cell型，单镜头各帧的wand像素坐标
%                                       每个cell: 1*2pointN    单镜头单帧的wand像素坐标，格式为[u1,v1,u2,v2,u3,v3,...]
% inexI     camN*10  各镜头的内外参数，每行为一个镜头的参数.[ f dx u0 v0 tx ty tz rx ry rz ]
% shield    1*camN   cell型，各镜头的屏蔽区域
%                    每个cell： n*4，单镜头屏蔽区域的斜对角点坐标，一行表示一块矩形区域，
%                    [右下角u, 右下角v, 左上角u, 左上角v]
% output：
% uv_useful frameN*1 cell型，各帧筛选后的wand坐标
%                    每个cell： n*5，单帧中筛选过后的wand坐标，每行为一个镜头，格式为[ 镜头号, u1, v1, u2, v2 ]
% FrameUsed camN*1   各镜头的可用帧数


camN = length(wand) ;     % 使用的镜头个数

frameN = zeros(1,camN) ;  %每个镜头的帧数
for icam = 1:length(wand)
    frameN(icam) = length(wand{icam}) ; 
end

useM = zeros(camN,11) ;   % 使用的镜头M矩阵
for i = 1:camN
    useM(i,:) = buildM(inexI(i,:)) ;
end

wand_Shield = cell(1,camN) ;         % 屏蔽后的wand坐标，按镜头存储
uv_useful = cell(min(frameN),1) ;    % 筛选后的wand坐标，按帧存储
FrameUsed = zeros(camN,1);           % 屏蔽筛选后每个镜头的可用帧数，以后可能会用的到……

for iframe = 1:min(frameN)
    uvtemp = zeros(1,camN*4) ; % 一帧中所有镜头点的uv坐标，排成一行的形式
    for icam = 1:camN 
        % 屏蔽杂点
        wand_Shield{icam}{iframe} = ShieldPoint(wand{icam}{iframe},shield{icam}) ;
        
        % 第一步筛选，去掉不等于3个点以及3个点不共线的图像
        uvtemp(icam*4-3:icam*4) = wand_screen_LM_1(wand_Shield{icam}{iframe}) ;
    end
    % 第二步筛选，去掉重建和重投影坐标偏离过大的图像
    uv_useful{iframe} = wand_screen_LM_2(uvtemp,useM) ;
    FrameUsed(uv_useful{iframe}(:,1)) = FrameUsed(uv_useful{iframe}(:,1)) + 1 ;
    
end

% 去掉数量不足2个镜头的帧
i = 1 ;
while i < length(uv_useful)
    if 2 > size(uv_useful{i},1)  
        uv_useful(i) = [] ; 
    else
        i = i + 1; 
    end
end

savepath = ['input\wand_shield'] ;
save(savepath,'wand')
fprintf('屏蔽杂点后的wand文件已保存在 %s \n', savepath)

save input\uv_useful uv_useful
disp('第二步筛选后的wand文件已保存在 input\uv_useful ')
