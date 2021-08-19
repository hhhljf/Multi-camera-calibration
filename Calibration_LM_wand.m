function [d, M0, uvi_useful ] = Calibration_LM_wand(inexIk, uv_useful)

% LM优化中求解的光束法平差目标函数（捆绑调整）
% 按帧计算，重建3D坐标，计算杆长误差；
% 用重建的3D坐标重投影到各摄像头，计算u,v方向上的误差.
% 含畸变，畸变方向：实际坐标 -(畸变)-> 理想坐标
% u,v 方向上的误差用理想坐标计算，整个流程为：
% 实际坐标 -(畸变)-> 理想坐标 -(重建)-> 3D坐标 -(投影)-> 理想坐标'
%                                  比较3D杆长误差      比较归一化平面上的2D误差

% input:
% inexIk        1*12camN    [inexI, k]，内外参加畸变系数 [ f dx u0 v0 tx ty tz rx ry rz k1 k2 ] 
% uv_useful     frameN*1    cell 每个cell中存储一帧中通过的镜头号及a,b两点的uv坐标，
%                           每行格式[ 镜头号, ua, va, ub, vb ]
% indparameter    整数      index parameter 计算Jacobian矩阵时发生改变的参数的下标
%                           -1 表示没有参数发生改变
% 
% output:
% d             1列          2D、3D误差向量
% M0            camN*11      本次计算后的M矩阵
% uvi_useful     frameN*1     畸变校正后的理想坐标

parameterN = 12 ; %每个镜头的参数个数
uvi_useful = cell(size(uv_useful)) ; %畸变校正后的理想坐标

% if size(inexIk,2)==1, inexIk=inexIk(:); end

frameN = length(uv_useful) ; %总帧数
camN = length(inexIk)/parameterN ;%camN 表示镜头的个数

inexIn = zeros(1,camN*10) ; % camN个镜头的内外参
kn = zeros(1,camN*2) ; %camN个镜头的畸变系数
for i = 1:camN
    inexIn(10*i-9:10*i) = inexIk(12*i-11:12*i-2) ;
    kn(2*i-1:2*i) = inexIk(12*i-1:12*i) ;
end

% 用inexIn生成M
M0 = zeros(camN,11) ;
for i = 1:camN
    M0(i,:) = buildM(inexIn(10*i-9:10*i)) ;
end

% time_dis = 0 ;          %畸变矫正耗时
% time_rebuild = 0 ;      %重建计时
% time_reprojection = 0 ; %重投影计时

% 按帧计算，重投影的像素uv差
d = []; %zeros( frameN*(1+), 1 ) ;
% d = zeros(frameN,1) ; %试试看用和呢……

for iframe = 1:frameN
    uv = uv_useful{iframe} ;
    M = M0( uv(:,1), :) ;
% tic    
    % 畸变校正
    uvi = zeros(size(uv)) ; % uv_ideal 去掉畸变的理想像素坐标
    uvi(:,1) = uv(:,1) ;
    for iuv = 1:size(uv,1) ;
        uvi(iuv,2:5) = adddistortion(uv(iuv,2:5), inexIn(uv(iuv,1)*10-9:uv(iuv,1)*10),  kn(uv(iuv,1)*2-1:uv(iuv,1)*2)) ; %加畸变
    end
    uvi_useful{iframe} = uvi ;
% time_dis = time_dis + toc ;    

%     usefulNum = size(uv,1) ;
%     % 两两重建，计算a,b点的3D坐标
%     k = 1 ;
%     xyza = zeros(usefulNum*(usefulNum-1)/2, 3) ; % a点对应的 [ x, y, z ]
%     xyzb = zeros(usefulNum*(usefulNum-1)/2, 3) ; % b点对应的 [ x, y, z ]
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
    % 两两重建实在是太慢了……直接用所有镜头最小二乘吧……
    tempuv = uvi(:,2:3)' ;
    tempuv = tempuv(:)' ;
    mxyza = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
    tempuv = uvi(:,4:5)' ;
    tempuv = tempuv(:)' ;
    mxyzb = rebulid_3D_UnspecCam_LM2(tempuv,M) ;
% %     d = [d; 38.7298*(500-norm(mxyza-mxyzb))] ;
    d = [d; 500-norm(mxyza-mxyzb)] ;
%     d(iframe) = abs(500-norm(mxyza-mxyzb)) ; %试试看用和呢……
% time_rebuild = time_rebuild + toc ;    

    % 重投影
% tic
    for i = 1:size(uv,1)
        uva = reprojection( M(i,:), mxyza ) ;
        uvb = reprojection( M(i,:), mxyzb ) ;
        d = [d; (uvi(i,2:5)-[uva,uvb])'   ] ;
%         d = [d; norm(uvi(i,2:3)-uva); norm(uvi(i,4:5)-uvb) ] ; % d49569
%         d(iframe) = d(iframe) + sum(abs((uvi(i,2:5)-[uva,uvb]))) ; % d2981
    end
    
%     d=sum(abs(d)); %试试看用和呢……
    
% time_reprojection = time_reprojection + toc ;
end

% fprintf('畸变校正耗时：%f \n',time_dis)
% fprintf('重建影耗时：%f \n',time_rebuild)
% fprintf('重投影耗时：%f \n',time_reprojection)

% ------------------------------------------------
% % 加入L框架的信息
% 
% xyz=[0	0	0 
% 200	0	0
% 600	0	0
% 0	400	0
% ]; % 规定好的L型框架世界坐标
% 
% Luv = load('C:\Users\Boat\Desktop\实采\data_uvc\CalFrame1_shield') ; %直接读取转换好的L坐标
% % Luv = load('input\Luv') ;
% Luv = Luv.uv ;
% iframe = randperm(length(Luv{1}),1) ;
% % iframe = randperm(min(frameN),1) ;  %随便选一帧。。。
% Luvr = zeros(camN,8) ; % Luv_real 各镜头实际拍到的L坐标
% for icam = 1:camN
%     Luvr(icam,:) = Luv{icam}{iframe} ;
% end
% 
% % 计算重建的3D误差
% Luv = zeros(4,2*camN) ; %各镜头的L坐标转换为标定函数接口的格式
% for i = 1:camN
%     for k = 1:4
%         Luv(k,2*(i-1)+1:2*i) = Luvr(i,2*(k-1)+1:2*k) ;
%     end
%     Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
% end
% 
% dd = zeros(4,1) ; % L型框架四个点的3D误差
% for k = 1:4
%     dd(k) = norm( xyz(k,:) - rebulid_3D_UnspecCam_LM(Luv(k,:),M0) ) ;
% end
% 
% d = [d; sqrt(iframe/15)*dd] ;

% 计算重投影的2D误差
% Luvi = zeros(camN,8) ; % Luv_ideal 各镜头重投影的L坐标
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


