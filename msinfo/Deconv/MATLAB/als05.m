%ALS(繰り返し計算の外でスペクトルの規格化−横長)
%Xを規格化
clear all

load nfhplc
X=[d1;d2;d3;d4];
% multiblock
mb=4;
bs=size(d1,1);

%load dadstart1


lambda1=eps;% 影響ないほう
%lambda2=0.00001;% 影響ないほう8.5700e+010
lambda2=eps;
% optimal?=0.00000000000001

%成分数
com=4;

for j=1:100% iter
   
C=rand(size(X,1),com);
C0=C;

for k=1:10000
   
   A=inv(C'*C+lambda1*eye(size(C',1)))*C'*X;
   %A=inv(C'*C)*C'*X;
   A=max(eps,A);
   
   % 規格化
   for i=1:com
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A'+lambda2*eye(size(A,1)));
   C=max(eps,C);
   
   
   % unimodality constraint using MCR-ALS toolbox
   rmod=0.1;
   C=unimod(C,rmod,1);

   
end

%plot(C,'k')
%hold on

L(j)=norm(X-C*A);

if norm(X-C*A)<10
   break
end


end


%------------------
%figure(1),
%for i=1:com
%   subplot(com,1,i),plot(s(i,:))
%end

%figure(2),
%for i=1:com
%   subplot(com,1,i),plot(A(i,:))
%end

%plot(C)

%plot(C(1:51,:))

norm(X-C*A)
plot(A')
j



