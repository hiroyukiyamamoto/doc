
clear all

for kurikaeshi=1:100
   
   
   
   
%---------------
% data
%---------------
%load XX.mat
%X=XX';








load Y


data=Y;
X=Y;

[r,c]=size(X);

%----------------------
% �����ݒ�l
%----------------------
maxIter=1000;% �J��Ԃ���

com=2;% ������
eta=0.1;% �w�K��

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
X=M;% �V�����f�[�^�Ƃ��čl����

%--------------------------------------------------------------
% whitening Z=VX;
Xm=(X'-ones(c,1)*mean(X'))';% ���ϒl��������X

%-----------------------
%Xm=X-ones(com,1)*mean(X);
%Xm(1,:)=mean(X);% ���ς��1�听��
%for i=2:com
%   Xm(i,:)=[Xm(i,:)];% �ȍ~����i�听��
%end
%------------------------
%----------------------------------------------------------------

%figure(3),plot(Xm(1,:),Xm(2,:),'*')

C=Xm*Xm';
V=C^(-1/2);
Z=V*X;

X=Z;
%----------

% �w�K�̔���
for i=1:maxIter
   y=W*X;
   yp=min(0,y);% f(y)
   W=W-eta*(yp*y'-y*yp')*W;
end

G=max(0,y');


T=Y*y';% score
F=T';








X=Y';

[m,n]=size(X);
com=2;

%G=rand(m,com);
%F=rand(com,n);

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
