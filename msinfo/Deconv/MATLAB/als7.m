%ALS(�J��Ԃ��v�Z�̊O�ŃX�y�N�g���̋K�i���|����)
%X���K�i��
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
   
   % �K�i��
   for i=1:com
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A');
   C=max(0,C);
   
   E(k)=norm(X-C*A);
end

plot(C)
