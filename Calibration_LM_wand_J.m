function d = Calibration_LM_wand_J(inexIk, uv_useful, uvi_useful, M0,  indparameter)

% LM优化中求解的光束法平差目标函数（捆绑调整）
% 计算Jacobian矩阵专用……
% 与 Calibration_LM_wand 不同的是，只计算参数发生改变的镜头的M矩阵和理想坐标，以提升效率.
% 
% 按帧计算，重建3D坐标，计算杆长误差；
% 重建的3D坐标重投影到各摄像头，计算u,v方向上的误差.
% 含畸变，畸变方向：实际坐标 -(畸变)-> 理想坐标
% u,v 方向上的误差用理想坐标计算，整个流程为：
% 实际坐标 -(畸变)-> 理想坐标 -(重建)-> 3D坐标 -(投影)-> 理想坐标'
%                                  比较3D杆长误差      比较归一化平面上的2D误差

% input:
% inexIk        1*12camN    [inexI, k]，内外参加畸变系数 [ f dx u0 v0 tx ty tz rx ry rz k1 k2 ] 
% uv_useful     frameN*1    cell 每个cell中存储一帧中通过的镜头号及a,b两点的uv坐标，
%                           每行格式[ 镜头号, ua, va, ub, vb ]
% uvi_useful    frameN*1    参数改变前，畸变校正过的理想像素坐标，与uv_useful对应
% indparameter    整数      index parameter 计算Jacobian矩阵时发生改变的参数的下标
%                           -1 表示没有参数发生改变
% M0            camN*11      参数改变前，各镜头的M矩阵，一行表示一个镜头的参数
% 
% output:
% d             1列          2D、3D误差向量

parameterN = 12 ; %每个镜头的参数个数

if size(inexIk,2)==1, inexIk=inexIk(:); end

% 计算发生改变的是哪个镜头的哪个参数
dpar  = rem(indparameter,parameterN) ; %发生改变的是第几个参数
if dpar==0
    dpar  = parameterN ;
    dcam = fix(indparameter/parameterN) ; %发生改变的镜头号
else
    dcam = fix(indparameter/parameterN)+1 ;
end


frameN = length(uv_useful) ; %总帧数
camN = length(inexIk)/12 ;%camN 表示镜头的个数

inexIn = zeros(1,camN*10) ; % camN个镜头的内外参
kn = zeros(1,camN*2) ; %camN个镜头的畸变系数
for i = 1:camN
    inexIn(10*i-9:10*i) = inexIk(12*i-11:12*i-2) ;
    kn(2*i-1:2*i) = inexIk(12*i-1:12*i) ;
end

% 若发生改变的是内外参，则需要重新计算该镜头的M矩阵
if dpar<11
    M0(dcam,:) = buildM(inexIn(10*dcam-9:10*dcam)) ; % 用inexIn生成M
end

% time_dis = 0 ;          %畸变矫正耗时
% time_rebuild = 0 ;      %重建计时
% time_reprojection = 0 ; %重投影计时

% 按帧计算，重投影的像素uv差
d = []; %zeros( frameN*(1+), 1 ) ;
% d = zeros(frameN,1) ; %试试看用和呢……
for iframe = 1:frameN

% tic    
    % 畸变校正
    uvi = uvi_useful{iframe} ; % uv_ideal 去掉畸变的理想像素坐标
    M = M0( uvi(:,1), :) ;
    
    if dpar>10 %如果是畸变参数发生改变，则需要重新计算该镜头的理想坐标
        temp = find(uvi(:,1)==dcam) ;
        if ~isempty(temp)
            uv = uv_useful{iframe} ; % 含畸变的实际像素坐标
            uvi(temp,2:5) = adddistortion(uv(temp,2:5), inexIn(uv(temp,1)*10-9:uv(temp,1)*10),  kn(uv(temp,1)*2-1:uv(temp,1)*2)) ; %加畸变
        end
    end
% time_dis = time_dis + toc ;    

% tic
    % 两两重建实在是太慢了……直接用所有镜头最小二乘吧……
    tempuv = uvi(:,2:3)' ;
    tempuv = tempuv(:)' ;
    mxyza = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    tempuv = uvi(:,4:5)' ;
    tempuv = tempuv(:)' ;
    mxyzb = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    d = [d; 500-norm(mxyza-mxyzb)] ;
%    d(iframe) = abs(500-norm(mxyza-mxyzb)) ; 
%     d(iframe*2-1) = 500-norm(mxyza-mxyzb) ;  %试试看用和呢……
% time_rebuild = time_rebuild + toc ;   

    % 重投影
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

% d=sum(abs(d)); %试试看用和呢……

% fprintf('畸变校正耗时：%f \n',time_dis)
% fprintf('重建影耗时：%f \n',time_rebuild)
% fprintf('重投影耗时：%f \n',time_reprojection)


end
