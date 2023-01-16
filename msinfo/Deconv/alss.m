clear all

load yosi
C=a;

[m,n]=size(X);


for t=1:100
   
   %A=pinv(C)*X;% X=CA
   A=mcdregres(C,X);
   A=A.fitted;
   A=max(0,A);
   
   %C=X*pinv(A);
   C=mcdregres(A,X);
   C=C.fitted;
   C=max(0,C);
   
   E(t)=norm(X-C*A);
   
   %if t>2 & abs(E(t)-E(t-1))<0.001
   %   break
   %end
   
   %if E(t)<0.1
   %   break
   %end
   
end

%min(E)
%figure(4),plot(x,A(1,:),'k')
%hold on;figure(4),plot(x,A(2,:),':k')

figure(1),plot(C)
figure(2),plot(A')
figure(3),plot(E)
min(E)
