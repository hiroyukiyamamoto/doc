clear all

load Y
%Y=Y';
for i=1:100
   i
   
   C=rand(size(Y,1),2);
   %C=[100*rand(size(Y,1),1) rand(size(Y,1),1)];

lambda=0;
%lambda=0;

[m,n]=size(Y);

for t=1:100
   %lambda=lambda0/t;
      
   A=inv(C'*C+lambda*eye(size(C,2)))*C'*Y;
   %M1(t)=cond(C'*C);
   A=max(0,A);
   
   C=Y*A'*inv(A*A'+lambda*eye(size(A,1)));
   %M2(t)=cond(A*A');
   C=max(0,C);

   %E(t)=norm(Y-C*A);
end

%a(i)=min(M2);

%figure(4),plot(C(:,1),'k')
%hold on;figure(4),plot(C(:,2),':k')
%figure(4),plot(A(1,:)/norm(A(1,:)),'k')
%hold on;figure(4),plot(A(2,:)/norm(A(2,:)),':k')
figure(4),plot(A(1,:),'k')
hold on;figure(4),plot(A(2,:),'k')


%figure(5),%subplot(3,1,1),plot(E,'k')
%subplot(2,1,1),plot(M1,'k')
%subplot(2,1,2),plot(M2,'k')
end
