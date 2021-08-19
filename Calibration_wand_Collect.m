function [ uv_useful, FrameUsed ] = Calibration_wand_Collect(wand, inexI, shield)

% wand�궨�е�ɸѡ����
% ��Ҫ���������ӵ㡢��һ��ɸѡ�͵ڶ���ɸѡ������������
% 
% input:
% wand      1*camN   cell�ͣ�����ͷ��wand�������꣬����cell�ṹ����
%                    ÿ��cell: 1*frameN cell�ͣ�����ͷ��֡��wand��������
%                                       ÿ��cell: 1*2pointN    ����ͷ��֡��wand�������꣬��ʽΪ[u1,v1,u2,v2,u3,v3,...]
% inexI     camN*10  ����ͷ�����������ÿ��Ϊһ����ͷ�Ĳ���.[ f dx u0 v0 tx ty tz rx ry rz ]
% shield    1*camN   cell�ͣ�����ͷ����������
%                    ÿ��cell�� n*4������ͷ���������б�Խǵ����꣬һ�б�ʾһ���������
%                    [���½�u, ���½�v, ���Ͻ�u, ���Ͻ�v]
% output��
% uv_useful frameN*1 cell�ͣ���֡ɸѡ���wand����
%                    ÿ��cell�� n*5����֡��ɸѡ�����wand���꣬ÿ��Ϊһ����ͷ����ʽΪ[ ��ͷ��, u1, v1, u2, v2 ]
% FrameUsed camN*1   ����ͷ�Ŀ���֡��


camN = length(wand) ;     % ʹ�õľ�ͷ����

frameN = zeros(1,camN) ;  %ÿ����ͷ��֡��
for icam = 1:length(wand)
    frameN(icam) = length(wand{icam}) ; 
end

useM = zeros(camN,11) ;   % ʹ�õľ�ͷM����
for i = 1:camN
    useM(i,:) = buildM(inexI(i,:)) ;
end

wand_Shield = cell(1,camN) ;         % ���κ��wand���꣬����ͷ�洢
uv_useful = cell(min(frameN),1) ;    % ɸѡ���wand���꣬��֡�洢
FrameUsed = zeros(camN,1);           % ����ɸѡ��ÿ����ͷ�Ŀ���֡�����Ժ���ܻ��õĵ�����

for iframe = 1:min(frameN)
    uvtemp = zeros(1,camN*4) ; % һ֡�����о�ͷ���uv���꣬�ų�һ�е���ʽ
    for icam = 1:camN 
        % �����ӵ�
        wand_Shield{icam}{iframe} = ShieldPoint(wand{icam}{iframe},shield{icam}) ;
        
        % ��һ��ɸѡ��ȥ��������3�����Լ�3���㲻���ߵ�ͼ��
        uvtemp(icam*4-3:icam*4) = wand_screen_LM_1(wand_Shield{icam}{iframe}) ;
    end
    % �ڶ���ɸѡ��ȥ���ؽ�����ͶӰ����ƫ������ͼ��
    uv_useful{iframe} = wand_screen_LM_2(uvtemp,useM) ;
    FrameUsed(uv_useful{iframe}(:,1)) = FrameUsed(uv_useful{iframe}(:,1)) + 1 ;
    
end

% ȥ����������2����ͷ��֡
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
fprintf('�����ӵ���wand�ļ��ѱ����� %s \n', savepath)

save input\uv_useful uv_useful
disp('�ڶ���ɸѡ���wand�ļ��ѱ����� input\uv_useful ')
