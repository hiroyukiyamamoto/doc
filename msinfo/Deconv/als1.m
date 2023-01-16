clear all

%C0=rand(5,2)
%A0=rand(2,5)
%X=C0*A0;
%X=X+0.05*randn(5,5);

load X

%C=rand(5,2);
%A=rand(2,5);
%X=C*A;

%load X
 
 %X=X+0.1*randn(5,5);

% å≥ÇÕ2ê¨ï™ÇÃCÇ∆SÇä|ÇØçáÇÌÇπÇΩÇ‡ÇÃ

% ê¨ï™êî
com=2;
lambda=0;

C=[0.6457 0.3349
   0.7863 0.3182
   0.4057 0.3690
   0.8828 0.9867
   0.5118 0.2563];

%C=[0.0430 0.9438
%   0.3655 0.6551
%   0.8067 0.6509
%   0.3627 0.7282
%   0.9077 0.2389];

%C=rand(5,2);
%load C

%load CC
%C=CC;


for k=1:1
   
   A=inv(C'*C)*C'*X;
   A=max(0,A);
   

      
   %C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
   %C=max(0,C);
   C(:,1)=X*A(:,1)'*inv(A(:,1)*A(:,1)'+lambda*eye(size(A(:,1),1)));
   C(:,2)=X*A(:,2)'*inv(A(:,2)*A(:,2)'+lambda*eye(size(A(:,2),1)));
   
  
   
   E(k)=norm(X-C*A,'fro');
   
end

A
C


