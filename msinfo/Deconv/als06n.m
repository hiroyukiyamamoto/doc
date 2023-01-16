%ALS(ŒJ‚è•Ô‚µŒvZ‚ÌŠO‚ÅƒXƒyƒNƒgƒ‹‚Ì‹KŠi‰»|‰¡’·)
%X‚ğ‹KŠi‰»
clear all

load Y
X=Y;


lambda=eps;% ‰e‹¿‚È‚¢‚Ù‚¤

%¬•ª”
com=2;

for h=1:100
   h

C0=rand(size(X,1),com);
  C=C0;

  for k=1:100
     
      A=pinv(C)*X;
      %A=inv(C'*C)*C'*X;
      A=max(eps,A);
      
            A=max(eps,A);
   rmod=0.1;A=unimod(A,rmod,1);
      
      %for i=1:com
      %   A(i,:)=A(i,:)/norm(A(i,:));
      %end
         
   C=X*pinv(A);
   %C=X*A'*inv(A*A');
   C=max(eps,C);
   

   
end

E(h)=norm(X-C*A);
plot(A','k')
hold on;
end

%subplot(2,1,1),plot(GA)
%subplot(2,1,2),plot(GC)


%figure(2),
%subplot(3,1,1),plot(A(1,:)','k')
%subplot(3,1,2),plot(A(2,:)','k')
%subplot(3,1,3),plot(A(3,:)','k')
%figure(2),
%subplot(2,1,1),plot(C(:,1))
%subplot(2,1,2),plot(C(:,2))
%E(end)



