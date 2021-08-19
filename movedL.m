function inexI = movedL(xf, Luv, L)

% ����������ϵƽ����ת����L���Ϊ��׼������ϵ��
% input��
% xf    10camN*1    LM�Ż�������о�ͷ�������
% Luv   4*2camN     L�Ϳ�ܵ�4���������о�ͷ�е��������꣬ÿ4*2������Ϊһ����ͷ��uv���ص�����
%                   ��������uv������˳����Բ���xyz�еĵ��ϸ��Ӧ�� 
% L     4*3         L�Ϳ������������ϵ�µ����꣬ͨ��Ϊ
%                   L=[   0	  0	0
%                       200	  0	0
%                       600	  0	0
%                         0	400	0 ];
% �������̷�Ϊ����
% 1.ƽ�ƣ�����ǰ��������ؽ���L�Ϳ�ܵ�������ʵ���趨��L�Ϳ������Ƚϣ�
%         ���ؽ�����L�Ϳ������ԭ��ƽ�Ƶ��趨��ԭ�㣬���о�ͷ��λ�Ʋ�����֮ƽ�Ƽ��ɣ�
% 2.��ת�����ڵ�ǰ����ͷ����ת��һ���Ƕȣ�ֱ����ԭ�ռ�����ת��ͷ���׽���ת��������Ϊ����ͷ��ŷ������ת��
%         ���¼�����ת������Ϊ��ֱ�ӣ���ת��������Ϊ������L�Ϳ����תexRd���趨λ�ã��ٰ����Ե���ת������תexR0��
%         ����������ת����exR=exR0*exRd���˴�exR0��exRd��λ�ò��ɽ���������������ת����exR������ת����.

% Luv = [1027	698	908.0999756	634	1023.400024	424	942.5	499.5	850.4000244	609	760	505	935.125	506.5	949	428.5	1118	587.5	852.5	498.5	932.555542	484	1007	528	995.375	495.5	1028.5	516	854.333313	662.5	955.833313	578.5
% 815.75	708.5	825.5	671	924.5	434.5	1046.214233	506	810.666687	631.5	845.5	538	959.5	542.5	935.7999878	470	945.416687	606.5	1110.666626	519.5	887.611084	556	916.5	547	969.7000122	510	960.1428833	527	991.083313	681.5	1003.75	605.5
% 943.5	730	870.7999878	690	935.625	452.5	934.083313	518.5	942.416687	634.5	891.2999878	555	972.5999756	562	811.7999878	489	1059.916626	614.5	972.833313	550.5	962.055542	570	941.4000244	564	1085.5	521.5	818.875	548.5	1062.099976	691	1029.5	619.5
% 898.666687	747.5	967.2000122	729	960.2000122	492	916.166687	561.5	725.583313	678.5	805.833313	595.5	863.5	575	928.416687	492.5	1028.5	629	900.833313	566.5	1114.400024	599.5	994.9000244	601	914.7999878	540	1081	590.5	1021.928589	767	942.666687	636.5
% ];
% L = [0 0 0;
%     200 0 0 ;
%     600 0 0 ;
%     0 400 0 ;
%     ] ;

% load input\xf
camN = length(xf)/10 ;
inexI = reshape(xf,10,camN) ;  % �� 1*10n ��������Ϊ 10*n �ľ���
inexI = inexI' ;

M = [];
for i = 1:camN
    temp = buildM(inexI(i,:)) ;
    M = [M temp'];
end
M0=M' ;

% ��Luv���� 0 1 2 3 �ŵ��˳������
for i = 1:camN
    Luv(:,2*(i-1)+1:2*i) = Lframe( Luv(:,2*(i-1)+1:2*i) ) ;
end

% �ؽ�����L�Ϳ�ܵ�����
pwc = zeros(4,3) ;
for k = 1:4
    pwc(k,:) = rebulid_3D_UnspecCam_LM2(Luv(k,:),M0)  ;
end
tran = pwc(1,:) ; %ƫ�ƾ���
pwc = pwc - repmat(pwc(1,:),4,1) ; 

% �÷��Գƾ�����޵������������ת����
A = zeros(3*4,3) ; 
for i = 1:4
    A(3*i-2:3*i,:) = [  0   L(i,3)+pwc(i,3)   L(i,2)+pwc(i,2) 
                     L(i,3)+pwc(i,3)      0   -L(i,1)-pwc(i,1)
                     L(i,2)+pwc(i,2)  L(i,1)+pwc(i,1)     0 ] ;
    B(3*i-2:3*i) = [L(i,1)-pwc(i,1); L(i,2)-pwc(i,2); L(i,3)+pwc(i,3)] ;
end
B=B' ;
temp = A\B ;
S = [0 -temp(3) -temp(2); temp(3) 0 -temp(1); temp(2) temp(1) 0];
exRd = (eye(3,3)+S)*inv((eye(3,3)-S)) ;  % ��ת���� R=(I+S)*(I-S)^(-1)

% ����λ�Ʋ���
inexI(:,5:7) = inexI(:,5:7) - repmat(tran,camN,1) ;
temp = exRd'*inexI(:,5:7)' ;
inexI(:,5:7) = temp' ;

% ������ת����
temp = inexI(:,8:10) ;
temp = temp*pi/180 ;
for i = 1:size(temp,1)
    rx = temp(i,1); ry = temp(i,2); rz = temp(i,3); 
    exR0=[cos(ry)*cos(rz) cos(ry)*sin(rz) -sin(ry)
    sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz) sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz) sin(rx)*cos(ry)
    cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz) cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz) cos(rx)*cos(ry)
    ];
    temp(i,:) = exR2rxryrz(exR0*exRd) ;  %���� exR0*exRd ��˳���ܱ䣬������L�Ϳ����תexRd���ٰ����Ե���ת������תexR0��˳��
end
inexI(:,8:10) = temp*180/pi ;

end

function r = exR2rxryrz(exR)
% ����ת����exR����ת˳�� rz -> ry -> rx���ֽ�Ϊrx,ry,rz
% input:
% exR   3*3     ��ת����˳�� rz -> ry -> rx
% output:
% r     1*3     [rx ry rz]

ry = -asin(exR(1,3)) ;
rx = acos(exR(3,3)/cos(ry));
% if ( rx_hat>0 ) rx_hat=-rx_hat ; end  % ��ͷ���Ǹ��ӣ�����ֻ����x����������ת. 

if sign(sin(rx)*cos(ry)) ~= sign(exR(2,3))
    ry = sign(ry)*pi - ry ; 
    rx = acos(exR(3,3)/cos(ry));
%     if ( rx_hat>0 ) rx_hat=-rx_hat ; end 
end

rz = asin(exR(1,2)/cos(ry));

if sign(cos(ry)*cos(rz)) ~= sign(exR(1,1))
    rz = sign(rz)*pi - rz ;
end

r = [rx ry rz ] ;

end

function orderuv = Lframe(uv) 
%�������ܣ�����Lframe�ϵ�4���㣬�����������г�Ϊ 0 1 2 3 �ŵ�˳�� .uv �� orderuv ����4*2��ʽ
% 0 1 2 ��x���ϣ�3 ��y���ϡ�0Ϊ����ԭ�㡣

% ����˼�룺������ȷ��һ��ֱ�ߡ�����������������ֱ�߾���ܽ�����ֵ����Ϊ��5�����Ե��ڣ���Ϊ��һ��ֱ���ϵġ�����ֱ�ߵ��ϵ���3�ŵ㡣
% ���� 0 1 2 ��֮��ľ��롣��������� 0 2 �㡣�ɴ�ȷ��������һ������1�ŵ㡣 ����1�ŵ���������0�ŵ㣬����Զ����2�ŵ㡣

% Ϊ�˱�����д�����е��ȷ����
uv0=uv(1,:) ; 
uv1=uv(2,:) ;
uv2=uv(3,:) ;
uv3=uv(4,:) ;

orderuv=zeros(4,2) ;
% d ��ʱ��������ʾ����

threshold =5 ;% ��ֵ�趨Ϊ5

%�� 0 1 ����ֱ�ߣ������ж�

d1=disPtoLine(uv0,[uv1 uv2]);
d2=disPtoLine(uv0,[uv1 uv3]);
d3=disPtoLine(uv0,[uv2 uv3]);

if d1>threshold && d2>threshold && d3>threshold
    orderuv(4,:)=uv0;
    orderuv(1:3,:)=Marker_sort([uv3;uv1;uv2]) ;
%     ��������
end


if d1<threshold
    orderuv(4,:)=uv3;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv2]) ;
%     ��������
end

    
    
if d2<threshold
    orderuv(4,:)=uv2;
    orderuv(1:3,:)=Marker_sort([uv0;uv1;uv3]) ;
%     ��������
    
end    
    
if d3<threshold
    orderuv(4,:)=uv1;
    orderuv(1:3,:)=Marker_sort([uv0;uv3;uv2]) ;
%     ��������
    
end

end %Lframe end



function New_marker=Marker_sort(marker)
%1,2,3�ŵ�����
New_marker = zeros(3,2);
d1_2=sqrt((marker(2,1)-marker(1,1))^2+(marker(2,2)-marker(1,2))^2); %1,2�ŵ�֮��ľ���
d2_3=sqrt((marker(2,1)-marker(3,1))^2+(marker(2,2)-marker(3,2))^2); %2,3�ŵ�֮��ľ���
d1_3=sqrt((marker(1,1)-marker(3,1))^2+(marker(1,2)-marker(3,2))^2); %1,3�ŵ�֮��ľ���
d = [ d2_3 d1_3 d1_2 ];
index = find(d==max(d));
New_marker(2,:) = marker( index,:) ; %��������� 0 2 �㡣�ɴ�ȷ��������һ������1�ŵ㡣
index = find(d<max(d)) ;
New_marker(1,:) = marker(find(d==max(d(index))),:) ; %����1�ŵ���������0�ŵ㣬
New_marker(3,:) = marker(find(d==min(d(index))),:) ; %����Զ����2�ŵ㡣
end %Marker_sort end

function [ d ] = disPtoLine( P, CorPts ) 
%�������ܣ�����㵽ֱ�ߵľ��루���������ά�ģ�
% ���룺corpts ��u1 v1 u2 v2����ʽ������ȷ��һ��ֱ�ߡ�
% P ��u v����ʽ
% �����d ��ʾ�㵽ֱ�ߵľ���

% Detailed explanation goes here 
if length(P) == 2
    l = [ CorPts(1) - CorPts(3), CorPts(2) - CorPts(4), 0 ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), 0 ]; 
elseif length(P) == 3
    l = [ CorPts(1) - CorPts(4), CorPts(2) - CorPts(5), CorPts(3) - CorPts(6) ]; 
    pl = [ P(1) - CorPts(1), P(2) - CorPts(2), P(3) - CorPts(3) ]; 
else
    error('������������~') ;
end
tem = cross(pl, l); 
d = norm( tem ) / norm( l ); 

end %disPtoLine end

