%ALS(åJÇËï‘ÇµåvéZÇÃäOÇ≈ÉXÉyÉNÉgÉãÇÃãKäiâªÅ|â°í∑)
%XÇãKäiâª
clear all

load nfhplc.mat
%X=X+abs(min(min(X)));
X=[d1;d2;d3;d4];

%ê¨ï™êî
com=4;
lambda=0;

%C=rand(size(X,1),com);
%C1=C;
load C1
C=C1;

for ite=1:1
   ite

   for k=1:500
      
      A=inv(C'*C)*C'*X;
      A=max(0,A);
      
      for i=1:com
         A(i,:)=A(i,:)/norm(A(i,:));
      end
            
      C=X*A'*inv(A*A'+lambda*eye(size(A*A')));
      C=max(0,C);
      
      C(1:51,:)=unimod(C(1:51,:),1,1);
      C(52:102,:)=unimod(C(52:102,:),1,1);
      C(103:153,:)=unimod(C(103:153,:),1,1);
      C(154:204,:)=unimod(C(154:204,:),1,1);
      
      E(k)=norm(X-C*A,'fro');
      NA(k)=norm(A,'fro');
      NC(k)=norm(C,'fro');
   
end

      plot(A','k')
      hold on
      
      N(ite)=E(end);
   end
   
   
  
   
   



%figure(2),plot(As,'k')
   

   