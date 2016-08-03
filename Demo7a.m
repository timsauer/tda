clear;len=5000;
t=(1:len);
n1=9;n2=4;a1=rand(n1,1);a2=rand(n2,1);
% X=[sin(.2*t+a(1));sin(.2*t+a(2));sin(.2*t+a(3));
%     1*sin(.2*sqrt(2)/3*t+b(1));1*sin(.2*sqrt(2)/3*t+b(2));1*sin(.2*sqrt(2)/3*t+b(3))];
f1=1;f2=sqrt(2);f2=1.28;
X=[sin(repmat(f1*t,n1,1)+repmat(a1,1,len));sin(repmat(f2*t,n2,1)+repmat(a2,1,len))];
% X=[sin(f1*t+a(1));sin(f1*t+a(2));sin(f1*t+a(3));
% sin(f2*t+b(1));sin(f2*t+b(2));sin(f2*t+b(3))];
X=X+.1*randn(size(X));
%Xpca2=3*randn(3,4*n)*X2;

%%% PCA Analysis %%%

 %[Xpca2,s] = PCA(X,6);
% 
% figure(11);
% plot3(Xpca(1,:),Xpca(2,:),Xpca(3,:));
% title('PCA (no delays)','Fontsize',18); 

%%% Delays + PCA %%%
X2 = [X(:,9:end); X(:,5:end-4); X(:,1:end-8)];
X2=dec(X2,300);
%[Xpca,s] = PCA(X2,10);
[Xpca,s] = PCA(X2,2*(n1+n2));
% 
figure(12);
plot3(Xpca(1,:),Xpca(2,:),Xpca(3,:)+Xpca(4,:),'.','MarkerSize',20);pause
title('PCA (after 3 delays)','Fontsize',18);
Q=[1 1 0 0;1 -1 0 0;0 0 1 1;0 0 1 -1];Q=Q/sqrt(2);
%Xpca(3:6,:)=Q*Xpca(3:6,:);
%[L,diam]=MPlex(Xpca([1 2 3 4 5 6],:),8);
[L,diam]=MPlex(Xpca,8);
