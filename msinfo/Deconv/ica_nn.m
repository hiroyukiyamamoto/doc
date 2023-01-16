% huu.m�����s
function [X,W]=ica_nn(x,eta,maxIter,com)
% non-negative ICA
% Algorithms for non-negative independent component analysis
% mark plumbley

% data matrix
%load Y
%X=Y;
X=x;

% �����ݒ�l

%maxIter=100;
% ��������3�̎��͖������Ȃ̂ŁA���s���邩������Ȃ�
%com=2;
%eta=0.1;

[m,n]=size(X);

%---------------------
% main
%---------------------

%---------------------- step1
W=eye(com,com);

% PCA�ɂ��f�[�^X��white��
[U,S,V]=svd(X,0);
M=V(:,1:com)';
% M:white�ɂ��ꂽ�f�[�^�s��X
% ���ϒl�������Ȃ�PCA��loading
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

X=X/norm(X);% �K�i��
% �J��Ԃ��r���ŋK�i������Ƃ��������Ȃ�
% �K�i������Ƃ��ɁA�S�����Ɏ����Ă����̂͑ʖڂ��H

plot(X');