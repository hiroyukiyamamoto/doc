clear all
%load oono
%Y=newsp108(1:3:108,:);
%Y=Y+1;

load Y
x=0.1:0.1:30;

[m,n]=size(Y);

A=rand(2,n);

for t=1:10
   
   C=Y*pinv(A);
   C=max(0,C);
   
   A=pinv(C)*Y;% Y=AC A=pinv(C)*Y
   A=max(0,A);
   
   E(t)=norm(Y-C*A);
end

figure(4),plot(C(:,1),'k')
hold on;figure(4),plot(C(:,2),':k')