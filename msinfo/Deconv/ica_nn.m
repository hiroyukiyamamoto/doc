% huu.mを実行
function [X,W]=ica_nn(x,eta,maxIter,com)
% non-negative ICA
% Algorithms for non-negative independent component analysis
% mark plumbley

% data matrix
%load Y
%X=Y;
X=x;

% 初期設定値

%maxIter=100;
% 成分数が3の時は未検討なので、失敗するかもしれない
%com=2;
%eta=0.1;

[m,n]=size(X);

%---------------------
% main
%---------------------

%---------------------- step1
W=eye(com,com);

% PCAによるデータXのwhite化
[U,S,V]=svd(X,0);
M=V(:,1:com)';
% M:whiteにされたデータ行列X
% 平均値を引かないPCAのloading
X=M;

for h=1:maxIter
   %-------------------------- step2
   Y=X;
   Yp=max(Y,0);
   Yn=min(Y,0);
   
   %-------------------------- step 3
   W=W+expm(-eta*(Yn*Yp'-Yp*Yn'));
   
   %--------------------------- step 8
   X=W*M;
end

X=X/norm(X);% 規格化
% 繰り返し途中で規格化するとおかしくなる
% 規格化するときに、全部正に持っていくのは駄目か？

plot(X');