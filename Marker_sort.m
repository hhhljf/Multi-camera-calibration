function New_marker=Marker_sort(marker)
%1,2,3�ŵ�����
New_marker = zeros(3,2);
d1_2=sqrt((marker(2,1)-marker(1,1))^2+(marker(2,2)-marker(1,2))^2); %1,2�ŵ�֮��ľ���
d2_3=sqrt((marker(2,1)-marker(3,1))^2+(marker(2,2)-marker(3,2))^2); %2,3�ŵ�֮��ľ���
d1_3=sqrt((marker(1,1)-marker(3,1))^2+(marker(1,2)-marker(3,2))^2); %1,3�ŵ�֮��ľ���
d = [ d2_3 d1_3 d1_2 ];
index = find(d==max(d));
if length(index)>1, New_marker=[]; return; end
New_marker(2,:) = marker( index,:) ; %��������� 0 2 �㡣�ɴ�ȷ��������һ������1�ŵ㡣
index = find(d<max(d)) ;

findmax = find(d==max(d(index))) ;
findmin = find(d==min(d(index))) ;
if length(findmax)==2, New_marker=nan(3,2); return; end

New_marker(1,:) = marker(find(d==max(d(index))),:) ; %����1�ŵ���������0�ŵ㣬
New_marker(3,:) = marker(find(d==min(d(index))),:) ; %����Զ����2�ŵ㡣
end %Marker_sort end
