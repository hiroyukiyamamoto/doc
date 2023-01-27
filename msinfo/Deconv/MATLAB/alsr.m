%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load dadstart2
%X=max(0,X);
X=X+abs(min(min(X)));


%¬•ª”
com=2;
lambda=0.0100;


%C=rand(size(X,1),com);
%load C
%load Y1
%C=Y1';

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
      
      E(k)=norm(X-C*A,'fro');
      NA(k)=norm(A,'fro');
      NC(k)=norm(C,'fro');
   
     
end

      plot(A','k')
      hold on
  