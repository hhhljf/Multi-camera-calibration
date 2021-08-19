function tuv = ShieldPoint(uv,shield)

% ���ε�������ͷ��֡�ĵ���uv�еĵ�
% input:
% uv        1*pointN    ������ͷ��֡�ĵ���uv���꣬��ʽΪ[ u1 v1 u2 v2 ... ]
% shield    rangeN*4    ���ε�������Ϣ��ÿ�б�ʾһ�����������б�Խǵ����꣬[���½�u, ���½�v, ���Ͻ�u, ���Ͻ�v]
% output:
% tuv       1*pointN'   ���κ�þ�ͷ��ʣ�µĵ����꣬��ʽΪ[ u1 v1 u2 v2 ... ]

tuv = [] ;
ituv = 1 ;
for i = 1:length(uv)/2
    flag = 0 ; %0��ʾ�����������ڣ�1��ʾ��
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




