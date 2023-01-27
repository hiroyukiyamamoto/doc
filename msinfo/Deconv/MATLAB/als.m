% Alternating least squares
% synthesis data

clear all

% data
load dadstart2
X=max(0,X);

% number of component
com=2;
lambda=0.0000001;
maxiter=100;

% initial value
C=rand(size(X,1),com);
%load C C


for k=1:maxiter
   
   A=inv(C'*C+lambda*eye(size(C'*C)))*C'*X;
   A=max(0,A);
   
   % normalized constraint
   for i=1:com
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A'+lambda*eye(size(A*A')));
   C=max(0,C);
   
   % Frobenious norm ||X-CA||
   E(k)=norm(X-C*A,'fro');
        
end

%plot(A','k')
%figure(2),plot(E,'k')