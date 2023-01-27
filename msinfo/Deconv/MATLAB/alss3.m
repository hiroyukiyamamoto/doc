clear all
%load oono
%Y=newsp108(1:3:108,:);
%Y=Y+1;
load Y
x=0.1:0.1:30;

[m,n]=size(Y);

C=rand(2,n);

for t=1:10
   
   A=Y*pinv(C);
   A=max(0,A);
   
   C=pinv(A)*Y;% Y=AC A=pinv(C)*Y
   C=max(0,C);
   
   E(t)=norm(Y-A*C);
end

figure(4),plot(x,C(1,:),'k')
hold on;figure(4),plot(x,C(2,:),':k')