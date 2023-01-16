%ALS(繰り返し計算の外でスペクトルの規格化－横長)
%Xを規格化
clear all

load dadstart3


X=max(X,0);% 0で切ると、負をそのままにした時とほとんど一緒の感じ
%X=X+abs(min(min(X)));% 正にすると、かなり変わる。ある意味では予想通り。
%X=Y;
%X=X+0.1;

%load Y
%X=Y;

%X=X(:,60:180);

%成分数
com=2;
lambda=0;


   
   %C=rand(size(X,1),com);
   load C

   
   
   for k=1:1
      
      A=inv(C'*C)*C'*X;
      A=max(0,A);
      
      %for i=1:com
      %   A(i,:)=A(i,:)/norm(A(i,:));
      %end
           
      C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
      inv(A*A'+lambda*eye(size(A,1)))
      C=max(0,C);

E(k)=norm(X-C*A,'fro');
end
   


figure(1),plot(A')
figure(2),plot(C)
  
  
   E2=norm(X-C*A,'fro');
   
 
   



%subplot(3,1,1),plot(E)
%subplot(3,1,2),plot(NC)
%subplot(3,1,3),plot(NA)
   