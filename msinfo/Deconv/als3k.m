%ALS(�J��Ԃ��v�Z�̊O�ŃX�y�N�g���̋K�i���|����)
%X���K�i��
clear all

load dadstart3


X=max(X,0);% 0�Ő؂�ƁA�������̂܂܂ɂ������ƂقƂ�ǈꏏ�̊���
%X=X+abs(min(min(X)));% ���ɂ���ƁA���Ȃ�ς��B����Ӗ��ł͗\�z�ʂ�B
%X=Y;
%X=X+0.1;

%load Y
%X=Y;

%X=X(:,60:180);

%������
com=2;
lambda=0;


   
   %C=rand(size(X,1),com);
   load C

   
   
   for k=1:1
      
      A=inv(C'*C)*C'*X;
      A=max(0,A);
      
      %for i=1:com
      %   A(i,:)=A(i,:)/norm(A(i,:));
      %end
           
      C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
      inv(A*A'+lambda*eye(size(A,1)))
      C=max(0,C);

E(k)=norm(X-C*A,'fro');
end
   


figure(1),plot(A')
figure(2),plot(C)
  
  
   E2=norm(X-C*A,'fro');
   
 
   



%subplot(3,1,1),plot(E)
%subplot(3,1,2),plot(NC)
%subplot(3,1,3),plot(NA)
   