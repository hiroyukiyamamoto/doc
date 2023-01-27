
clear all
load XX.mat

X=XX;

[n,m]=size(X);
%X=X+10*rand(r,s);

% scaling
%for j=1:m
%   XX(:,j)=X(:,j)/sqrt(var(X(:,j)));
%end

% kikakuka
%for j=1:m
%   XX(:,j)=X(:,j)/sum(X(:,j));
%end
%X=XX;

%------------

%--------------
% Šî’ê‚Ì”‚ğw’è
num=2;
%---------------------
r=num;

%--------------
%H=rand(num,m);% ‰Šú’l
%N=unprincomp(X');
%W=N(:,1:3);
W=rand(n,num);
W=W';

% kikakuka
for j=1:n
   WW(:,j)=W(:,j)/sum(W(:,j));
end
W=WW';

% scaling
%for j=1:num
%   W(:,j)=W(:,j)/sqrt(var(W(:,j)));
%end

V=X;
%H=rand(num,s);

for i=1:m
   H(:,i)=nnls(W,V(:,i));
end

% scaling
%for j=1:num
%   WWF(:,j)=W(:,j)/sum(W(:,j));
%end
%W=WWF;

%if itn='null'
%   itn=1;
%end

for iter=1:10
   
   WH=W*H;
   
   %for t=1:n
   %   for y=1:m
   %      if WH(t,y)<0
   %         WH=0;
   %      end
   %   end
   %end
   
   for i=1:n
      for a=1:r
         for u=1:m
            W1(i,a)=W(i,a)*sum(V(i,u)*H(a,u)/WH(i,u));%eq1
         end
      end
   end
   
   for i=1:n
      W2(i,:)=W1(i,:)/sum(W1(i,:));%eq2
   end
   
   W22(:,:,iter)=W2;
   
   %WH=W*H;
   
   for i=1:n
      for a=1:r
         for u=1:m
            H1(a,u)=H(a,u)*sum(W(i,a)*V(i,u)/WH(i,u));%eq3
         end
      end
   end

%HH=H1';
%for j=1:m
%   H2(j,:)=HH(j,:)/sum(HH(j,:));
%end
%H=H2';

W=W2;
H=H1;

%if rank(W')<num
%   break
%end

if isnan(W(1,1)) == 1
   break
end

%NN(iter)=norm(XX-W*H);
iter
end
bar(H')