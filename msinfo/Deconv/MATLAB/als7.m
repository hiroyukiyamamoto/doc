%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load ic2a
load dadstart3
X=X+abs(min(min(X)));
com=2;

%C=Yp';
C=ICZ';

for k=1:3
   
   A=inv(C'*C)*C'*X;
   A=max(0,A);
   
   % ‹KŠi‰»
   for i=1:com
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A');
   C=max(0,C);
   
   E(k)=norm(X-C*A);
end

plot(C)
