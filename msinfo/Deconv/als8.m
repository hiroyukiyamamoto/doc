%ALS(�J��Ԃ��v�Z�̊O�ŃX�y�N�g���̋K�i���|����)
%X���K�i��
clear all

%������
com=2;

C=[0.1 0.2
   0.2 0.8
   0.7 0.6
   0.5 0.4
   0.2 0.1
];

A=rand(2,100);
% �K�i��
for i=1:com
   A(i,:)=A(i,:)/norm(A(i,:));
end

X=C*A;
%load Yp
C=rand(5,2);

for k=1:100
   
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
