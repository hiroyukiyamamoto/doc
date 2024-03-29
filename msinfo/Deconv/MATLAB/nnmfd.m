

%----------------------------------
% non-negative matrix factorization
%
% m:サンプル数
% n:波長の数
% r:成分数
%
% W:基底（純スペクトル）
% H:展開係数（濃度行列）
%
% ---------------------------------

clear all


%--------------------
% data
%--------------------
load Y
V=Y';

[r,c]=size(V);
for i=1:r
   for j=1:c
      if V(i,j)==0
         V(i,j)=0.01;
      end
   end
end

%load XX.mat
%V=XX';
[n,m]=size(V);


%--------------------
% 基底の数を指定
%--------------------
num=3;
r=num;


%--------------------
% 初期値
%--------------------
W=rand(n,num);% 初期値
H=rand(num,m);


%-------------------
% main
%-------------------
for iter=1:1000
   
   
   %-------------------
   % 更新
   %-------------------
   WH=W*H;
   
   
   %----------------------
   % eq3
   %----------------------
   for i=1:n
      for a=1:r
         for u=1:m
            H1(a,u)=H(a,u)*sum(W(i,a)*V(i,u)/WH(i,u));%eq3
         end
      end
   end
   
   
   %---------------------
   % eq(sub)　
   % 文献には無いけれども、Hも規格化しておきたい！
   %---------------------
   for i=1:m
      for a=1:r
         H2(a,i)=H1(a,i)/sum(H1(:,i));%eq2
      end
   end
   
   
   %----------------------
   % 更新
   %----------------------
   H=H2;
   
   
   %--------------------
   % eq1
   %--------------------
   for i=1:n
      for a=1:r
         for u=1:m
            W1(i,a)=W(i,a)*sum(V(i,u)*H(a,u)/WH(i,u));%eq1
         end
      end
   end
   
   
   %---------------------
   % eq2
   %---------------------
   for i=1:n
      for a=1:r
         W2(i,a)=W1(i,a)/sum(W1(:,a));%eq2
      end
   end
   
   
   %-----------------------
   % 更新
   %-----------------------
   W=W2;
   
   %-----------------------
   % 更新
   %-----------------------
   %W=W1;
   
   
end


%---------------------
% 図示
%---------------------
%bar(W)

plot(W)

