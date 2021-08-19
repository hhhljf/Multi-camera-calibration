 function w = rebulid_3D_UnspecCam_LM2(uv,M)
 %����LM�Ż��е��ؽ���ȥ���˲в�ļ���������ٶȡ�
 %%%���룺uv,M
 % uv    Ϊ����Ķ����ͷ��2D������
 % ��ʽ��ÿ����Ϊһ����ͷ�����ݣ�m����ͷ����2m��  
 % M    Ϊ����Ķ����ͷ��M����
 % ��ʽ��ÿһ��Ϊһ����ͷ��M����n����ͷ��Ϊn*11�ľ���
  
 %%%�����w,s
 % w      Ϊ�ؽ���3D����
 %��ʽ��  uv�����һ���ؽ�һ����
 % s      Ϊ�����ͷ�����3D�в�
 %��ʽ��  uv�����һ�м���һ��ֵ
 
 uv_l = size(uv,1);%��ȡuv������
 M_l=size(M,1);   %��ȡM������
uv = uv' ;
 
 w=zeros(uv_l,3);
 for s_sol=1:uv_l
     
     %%����3D����
     A=zeros(2*M_l,3); b=zeros(2*M_l,1);%������A��b

% ֮ǰ�İ汾�����ܺ������ for icam = 1:M_l ������ͬ
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
 
    