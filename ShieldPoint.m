function tuv = ShieldPoint(uv,shield)

% 屏蔽掉单个镜头单帧拍到的uv中的点
% input:
% uv        1*pointN    单个镜头单帧拍到的uv坐标，格式为[ u1 v1 u2 v2 ... ]
% shield    rangeN*4    屏蔽的区域信息，每行表示一个矩形区域的斜对角点坐标，[右下角u, 右下角v, 左上角u, 左上角v]
% output:
% tuv       1*pointN'   屏蔽后该镜头中剩下的点坐标，格式为[ u1 v1 u2 v2 ... ]

tuv = [] ;
ituv = 1 ;
for i = 1:length(uv)/2
    flag = 0 ; %0表示不在屏蔽区内，1表示在
    for is = 1:size(shield,1) % index_shield
        if uv(i*2-1)<shield(is,1) & uv(i*2-1)>shield(is,3) & uv(i*2)<shield(is,2) & uv(i*2)>shield(is,4)
            flag = 1 ;
            break ;
        end
    end
    if flag==0 
        tuv(ituv*2-1:ituv*2) =  uv(i*2-1:i*2) ; 
        ituv = ituv + 1 ;
    end
end




