function Calibration_main()

% L框架 + wand 标定

% clear

% ======================================================================
% ==== Step1 L框架标定 ====
% 读入：手册内参，设定好的L型框架上4个点的世界坐标，L型框架上4个点在各镜头中的像素坐标
% 返回：各镜头的内、外、畸变参数初值，保存在input文件夹中
% （已在 Calibration_Lframe.m 中实现）

Calibration_Lframe()

% ======================================================================



% ======================================================================
% ==== Step2 wand优化 ====
% 读入：内、外、畸变参数初值，wand在各镜头中的像素坐标
% 返回：优化后的内、外、畸变参数、M矩阵，保存在input文件夹中
% （已在 Calibration_wand.m 中实现）

Calibration_wand()

% ======================================================================





