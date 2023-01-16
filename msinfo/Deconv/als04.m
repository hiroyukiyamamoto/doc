%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load nfhplc
X=[d1;d2;d3;d4];
mb=4;
bs=size(d1,1);


lambda=eps;% ‰e‹¿‚È‚¢‚Ù‚¤

%¬•ª”
com=4;

for j=1:52
   
   lambda=lambda*2;
   L(j)=lambda;
   
   C=rand(size(X,1),com);
   C0=C;
   
   for k=1:100
      
      A=inv(C'*C+lambda*eye(size(C',1)))*C'*X;
      A=max(eps,A);
      
      % ‹KŠi‰»
      for i=1:com
         A(i,:)=A(i,:)/norm(A(i,:));
      end
      
      C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
      C=max(eps,C);
      
      % unimodality constraint using MCR-ALS toolbox
      rmod=0.1;
      C(1:bs,:)=unimod(C(1:bs,:),rmod,1);
      C(bs+1:2*bs,:)=unimod(C(bs+1:2*bs,:),rmod,1); 
      C(2*bs+1:3*bs,:)=unimod(C(2*bs+1:3*bs,:),rmod,1);
      C(3*bs+1:4*bs,:)=unimod(C(3*bs+1:4*bs,:),rmod,1);
      
   end
   
   E(j)=norm(X-C*A);
   
end

plot(E,'o')

[a,b]=min(E);
L(b)


