clear all

load dadstart1

C=rand(size(X,1),2);
 
% å≥ÇÕ2ê¨ï™ÇÃCÇ∆SÇä|ÇØçáÇÌÇπÇΩÇ‡ÇÃ

% ê¨ï™êî
com=2;
lambda=0;


%C=rand(5,2);

for k=1:1000
   k
   
   A=inv(C'*C+lambda*eye(size(C,2)))*C'*X;
   A=max(eps,A);
   
   %for i=1:2
   %   A(i,:)=A(i,:)/norm(A(i,:));
   %end
   
   
   C'*C
   NC1(k)=norm(C(:,1));
   NC2(k)=norm(C(:,2));
   
   NC(k)=norm(C(:,2))/norm(C(:,1));
   
   C=X*A'*inv(A*A'+lambda*eye(size(A,1)));
   C=max(eps,C);
   
   E(k)=norm(X-C*A);
   
end

plot(E,'k')
