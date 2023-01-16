
function [W,H]=nnmf(Y,num)

X=Y';

[n,m]=size(X);
%X=X+10*rand(r,s);

% scaling
%for j=1:s
%   X(:,j)=X(:,j)/sqrt(var(X(:,j)));
%end

% kikakuka
%for j=1:s
%   X(:,j)=X(:,j)/sum(X(:,j));
%end

%------------

%--------------
% Šî’ê‚Ì”‚ğw’è
num=2;
r=num;
%---------------------

%--------------
W=rand(n,num);% ‰Šú’l
%N=unprincomp(X');
%W=N(:,1:3);

% kikakuka
%for j=1:num
%   W(:,j)=W(:,j)/sum(W(:,j));
%end

% scaling
%for j=1:num
%   W(:,j)=W(:,j)/sqrt(var(W(:,j)));
%end

V=X;
%H=rand(num,s);
H=pinv(W)*V;

%if itn='null'
%   itn=1;
%end

for iter=1:1
for i=1:n
   for a=1:r
      
      WH=W*H;
      
      for t=1:n
         for y=1:m
            if WH(t,y)<0.001
               WH(t,y)=0.001;
            end
         end
      end
      
      for u=1:r
         W(i,a)=W(i,a)*(sum((V(i,u)/WH(i,u))*H(a,u)));%eq1
      end
      
      for j=1:n
         W(i,a)=W(i,a)/sum(W(j,a));%eq2
      end
      
      for t=1:n
         for y=1:r
            if W(t,y)<0.001
               W(t,y)=0.001;
            end
         end
      end
      %for j=1:num
      %   W(:,j)=W(:,j)/sqrt(var(W(:,j)));
      %end
      
      WH=W*H;
      
      for t=1:n
         for y=1:m
            if WH(t,y)<0.001
               WH(t,y)=0.001;
            end
         end
      end
      
      for u=1:r
         H(a,u)=H(a,u)*sum(W(i,a)*(V(i,u)/WH(i,u)));%eq3
      end
            
   end

end

%HH=H';
%for j=1:m
%   H1(j,:)=HH(j,:)/sum(HH(j,:));
%end
%H=H1';

iter
end
plot(W)