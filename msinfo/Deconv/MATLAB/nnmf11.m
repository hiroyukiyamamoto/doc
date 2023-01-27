
clear all

for kurikaeshi=1:100
   
%---------------
% data
%---------------
%load XX.mat
%X=XX';
load Y

X=Y';

[m,n]=size(X);
com=2;

G=rand(m,com);
F=rand(com,n);

for iter=1:100
   
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
end


kurikaeshi
hold on;plot(G)

end
