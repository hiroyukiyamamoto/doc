%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load nfhplc.mat
%X=X+abs(min(min(X)));
X=[d1;d2;d3;d4];

lambda=10^-011;% ‰e‹¿‚È‚¢‚Ù‚¤
%¬•ª”
com=4;

for l=1:1
   l
   lambda=10^-011;
   
   for j=1:12
      j;
      
   
   L(j)=lambda;
   
   %C=rand(size(X,1),com);
   load C1
   C=C1;

   
   for k=1:10000
      
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



