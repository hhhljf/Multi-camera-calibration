 function w = rebulid_3D_UnspecCam_LM2(uv,M)
 %用于LM优化中的重建，去掉了残差的计算以提高速度。
 %%%输入：uv,M
 % uv    为输入的多个镜头的2D点数据
 % 格式：每两列为一个镜头的数据，m个镜头，有2m列  
 % M    为输入的多个镜头的M矩阵
 % 格式：每一行为一个镜头的M矩阵，n个镜头，为n*11的矩阵
  
 %%%输出：w,s
 % w      为重建的3D坐标
 %格式：  uv输入的一行重建一个点
 % s      为多个镜头计算的3D残差
 %格式：  uv输入的一行计算一个值
 
 uv_l = size(uv,1);%求取uv的行数
 M_l=size(M,1);   %求取M的行数
uv = uv' ;
 
 w=zeros(uv_l,3);
 for s_sol=1:uv_l
     
     %%计算3D坐标
     A=zeros(2*M_l,3); b=zeros(2*M_l,1);%不定长A和b

% 之前的版本，功能和下面的 for icam = 1:M_l 部分相同
%         for kk=1:3
%             for jj=1:M_l
%                 A(2*jj-1,kk)=uv(s_sol,2*jj-1)*M(jj,8+kk)-M(jj,kk);
%                 A(2*jj,kk)=uv(s_sol,2*jj)*M(jj,8+kk)-M(jj,4+kk);
%             end
%         end
%         for jj=1:M_l
%             b(2*jj-1,1)=M(jj,4)-uv(s_sol,2*jj-1);    
%             b(2*jj,1)=M(jj,8)-uv(s_sol,2*jj);
%         end

        for icam = 1:M_l
            A(icam*2-1,:) = uv(icam*2-1,s_sol) * M(icam,9:11) - M(icam,1:3) ;
            A(icam*2,:)   = uv(icam*2,s_sol)   * M(icam,9:11) - M(icam,5:7)  ;
            b(icam*2-1)   = M(icam,4)-uv(icam*2-1,s_sol) ;
            b(icam*2)     = M(icam,8)-uv(icam*2,s_sol) ;
        end
        
        w(s_sol,:)=(A\b)';
        
        
 end

 end
 
    