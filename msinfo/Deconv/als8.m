%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

%¬•ª”
com=2;

C=[0.1 0.2
   0.2 0.8
   0.7 0.6
   0.5 0.4
   0.2 0.1
];

A=rand(2,100);
% ‹KŠi‰»
for i=1:com
   A(i,:)=A(i,:)/norm(A(i,:));
end

X=C*A;
%load Yp
C=rand(5,2);

for k=1:100
   
   A=inv(C'*C)*C'*X;
   A=max(0,A);
   
   % ‹KŠi‰»
   for i=1:com
      A(i,:)=A(i,:)/norm(A(i,:));
   end
   
   C=X*A'*inv(A*A');
   C=max(0,C);
   
   E(k)=norm(X-C*A);
end

plot(C)
