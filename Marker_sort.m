function New_marker=Marker_sort(marker)
%1,2,3号点排序
New_marker = zeros(3,2);
d1_2=sqrt((marker(2,1)-marker(1,1))^2+(marker(2,2)-marker(1,2))^2); %1,2号点之间的距离
d2_3=sqrt((marker(2,1)-marker(3,1))^2+(marker(2,2)-marker(3,2))^2); %2,3号点之间的距离
d1_3=sqrt((marker(1,1)-marker(3,1))^2+(marker(1,2)-marker(3,2))^2); %1,3号点之间的距离
d = [ d2_3 d1_3 d1_2 ];
index = find(d==max(d));
if length(index)>1, New_marker=[]; return; end
New_marker(2,:) = marker( index,:) ; %距离最长的是 0 2 点。由此确定，另外一个点是1号点。
index = find(d<max(d)) ;

findmax = find(d==max(d(index))) ;
findmin = find(d==min(d(index))) ;
if length(findmax)==2, New_marker=nan(3,2); return; end

New_marker(1,:) = marker(find(d==max(d(index))),:) ; %距离1号点距离近的是0号点，
New_marker(3,:) = marker(find(d==min(d(index))),:) ; %距离远的是2号点。
end %Marker_sort end
