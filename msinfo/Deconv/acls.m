% Alternating constrained least squares
% for synthesis data

clear all

% data
load Y
X=Y;

% number of component
com=2;
lambda=0;
maxiter=1000;

% initial value
C=rand(size(X,1),com);

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

plot(A','k')