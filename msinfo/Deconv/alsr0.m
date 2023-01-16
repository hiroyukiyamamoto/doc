%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load dadstart2
%X=max(0,X);
X=X+abs(min(min(X)));

lambda=10^-011;% ‰e‹¿‚È‚¢‚Ù‚¤
%¬•ª”
com=2;

for l=1:1
   l
   lambda=10^-011;
   
   for j=1:12
      j;
      
   
   L(j)=lambda;
   
   %C=rand(size(X,1),com);
   load A2
   A=A2;
   
   for k=1:1000
      
      C=X*A'*inv(A*A'+lambda*eye(size(A*A')));
      C=max(0,C);
      
      C=unimod(C,1,1);
      
      A=inv(C'*C+lambda*eye(size(C'*C)))*C'*X;
      A=max(0,A);
      
      for i=1:com
         A(i,:)=A(i,:)/norm(A(i,:));
      end
      

      
      
      
      %E(k)=norm(X-C*A,'fro');
      
     
   end
   
   E1(j)=norm(X-C*A,'fro');

   
   lambda=lambda*10;% lambda‚ÌXV
   

end

%[a,b]=min(E);
%EE(l)=a;
%LL(l)=L(b);

   % ‹KŠi‰»




figure(l),plot(L,E1,'k-o')
hold on

end



