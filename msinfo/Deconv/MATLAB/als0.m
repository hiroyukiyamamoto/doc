%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load dadstart3
%X=Y;

lambda=0.00001;
l=0.05;

[m,n]=size(X);
E=eye(n,n);% n(”g’·”)~n‚Ì’PˆÊs—ñ
% second order differential operator
D=diff(E,1);

%¬•ª”
com=2;

%for h=1:100
C0=rand(size(X,1),com);
  C=C0;

   for k=1:100
      A=inv(C'*C)*(C'*X);
      %A=inv(C'*C)*C'*X;
      A=max(0,A);
   
   for i=1:size(A,1)
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
   C=(sgolayfilt(C,3,21));
   %C=X*A'*inv(A*A');
   C=max(0,C);
   E(k)=norm(X-C*A);
   
%   GA(k)=norm(A(1,:))/norm(A(2,:));
%   GC(k)=norm(C(:,1))/norm(C(:,2));
end
%GA(h)=norm(A(1,:))/norm(A(2,:));

%GC(h)=norm(C(:,1))/norm(C(:,2));

%plot(A')
%hold on;
%end

%subplot(2,1,1),plot(GA)
%subplot(2,1,2),plot(GC)

plot(A')



