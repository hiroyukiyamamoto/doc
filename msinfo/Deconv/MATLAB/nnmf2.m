
clear all

%---------------
% data
%---------------
load XX.mat
X=XX';

[m,n]=size(X);
com=2;

G=rand(m,com);
F=rand(com,n);

for iter=1:200
   
   T1=G'*X;
   T2=G'*G*F;
   T3=X*F';
   T4=G*F*F';
   
   for a=1:com      
    
      for u=1:n
         F(a,u)=F(a,u)*T1(a,u)/T2(a,u);
      end
      
     
      T1=G'*X;
      T2=G'*G*F;
      T3=X*F';
      T4=G*F*F';
      
      for i=1:m
         G(i,a)=G(i,a)*T3(i,a)/T4(i,a);
      end
      
      N(iter)=sum(sum(X-G*F));
      
   end
   
   for j=1:n
      F1(:,j)=F(:,j)/sum(F(:,j));
   end
   F=F1;
   
   
end

x=[27 28 29 39 40 41 42 43 51 53 54 55 56 67 68 69 79 81 82 84];
subplot(2,1,1),bar(x,G(:,1),'LineWidth',0.1)
subplot(2,1,2),bar(x,G(:,2),'LineWidth',0.1)