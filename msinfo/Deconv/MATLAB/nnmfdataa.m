
clear all

%---　元のデータ ---%
x=0.1:0.1:30;

m01=15;s01=3;m02=20;s02=1;m03=10;s03=2;
y1=exp(-(x-m01).^2./(2*s01*s01))./((sqrt(2*pi))*s01);
y2=exp(-(x-m02).^2./(2*s02*s02))./((sqrt(2*pi))*s02);
y3=exp(-(x-m03).^2./(2*s03*s03))./((sqrt(2*pi))*s03);

Ab=[y1;y2;y3];

%-- 設定値
m=5;% サンプル数
%--

n=3;% 固定値 
CoR=rand(m,n);
%CoR=abs(randn(m,n));

for i=1:m
   for j=1:n
      Co(i,j)=(CoR(i,j)/sum(CoR(i,:)))*1;
   end
end

clear m n

Y=Co*Ab;

%-------------------------------------------------------
X=Y';

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
% 基底の数を指定
num=3;
%---------------------
r=num;

%--------------
H=rand(num,m);% 初期値
%N=unprincomp(X');
%W=N(:,1:3);

% kikakuka
%for j=1:m
%   H1(:,j)=H(:,j)/sum(H(:,j));
%end
%H=H1;

% scaling
%for j=1:num
%   W(:,j)=W(:,j)/sqrt(var(W(:,j)));
%end

V=X;
%H=rand(num,s);

for i=1:n
   W(:,i)=nnls(H',V(i,:)');
end
W=W';

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
      for a=1:r
         W2(i,a)=W1(i,a)/sum(W1(:,a));%eq2
      end
   end
   
   if isnan(W2(1,1)) == 1
      break
   end
      
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
iter;
end

if isnan(W(1,1)) == 1
   break
end


