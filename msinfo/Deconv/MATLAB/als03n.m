%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all


load dadstart2
X=max(0,X);
%X=Y;


lambda=0;% ‰e‹¿‚È‚¢‚Ù‚¤

%¬•ª”
com=2;

for h=1:100
   
   h
  

C0=rand(size(X,1),com);
C=C0;
%load C

  for k=1:100
     
   A=inv(C'*C)*C'*X;
   A=max(0,A);
      
      for i=1:com
         A(i,:)=A(i,:)/norm(A(i,:));
      end
         
   C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
   C=max(0,C);
   %C=unimod(C,1,1);

   
end

   E=norm(X-C*A,'fro');
     
%E(h)=norm(X-C*A);
K(h)=E;

plot(C,'k')
hold on;
end

%subplot(2,1,1),plot(GA)
%subplot(2,1,2),plot(GC)


%figure(2),
%subplot(3,1,1),plot(A(1,:)','k')
%subplot(3,1,2),plot(A(2,:)','k')
%subplot(3,1,3),plot(A(3,:)','k')
%figure(2),
%subplot(2,1,1),plot(C(:,1))
%subplot(2,1,2),plot(C(:,2))
%E(end)

%plot(E,'o')

