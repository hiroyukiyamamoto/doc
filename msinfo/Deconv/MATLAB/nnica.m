clear all

load dadstart2
X=X+abs(min(min(X)));
X=X';

[n,m]=size(X);

com=2;
W=eye(com);

[E,D]=eig(cov(X'));
Z=D(end-com+1:end,end-com+1:end)^(-1/2)*E(:,end-com+1:end)'*X;

Y=W*Z;
eta=0.001;

for t=1:1000
   
   Yp=max(Y,0);
   Ym=min(Y,0);
   
   W=expm(-eta*(Ym*Yp'-Yp*Ym'))*W;
   Y=W*Z;
   
   [n,p]=size(X);
   O=(Z-W'*Yp).^2;
   P=sum(O(:));
   Ennr(t)=1/n*p*P;
   
end
plot(Ennr)

S=X*pinv(Yp);

Yp(1,:)=Yp(1,:)*norm(S(:,1));
Yp(2,:)=Yp(2,:)*norm(S(:,2));

S(:,1)=S(:,1)/norm(S(:,1));
S(:,2)=S(:,2)/norm(S(:,2));



figure(2),plot(y,S)
figure(3),plot(x,Yp')

