% non-negative ICA
% Algorithms for non-negative independent component analysis
% mark plumbley
clear all

% data matrix
load Y
data=Y;
X=Y;

% 初期設定値
maxIter=100;
% 成分数が3の時は未検討なので、失敗するかもしれない
com=2;
eta=0.1;

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

N=princomp(X);
M=N(:,1:com)';

for h=1:maxIter
   %-------------------------- step2
   Y=X;
   Yp=max(Y,0);
   Yn=min(Y,0);
   
   H=Yp*Yn'-Yn*Yp';
   G=Yn*Yp'-Yp*Yn';
   [p,q]=size(G);
   
   for i=2:p
      for j=1:i-1
         G(i,j)=0;
      end
   end
   
   if norm(G)<0.001
      break
   end
   
   l=1;
   for tau=0:0.1:100
      R=expm(tau*H);
      Y=R*X;
      
      Yp=max(Y,0);
      Yn=min(Y,0);
      J(l)=(1/2)*trace(Yn'*Yn);
      
      l=l+1;
   end
   
   yoko=0:0.1:100;
   %figure(h),plot(yoko,J)
   
   [m,n]=min(J);
   tau=0+0.1*(n-1);
   
   R=expm(tau*H);
   
   W=R*W;
   X=R*M;
   X=X/norm(X);
end


plot(X');