clear;len=3000;dt=1;
t=(1:len);
n1=5;n2=3;X=[];
for i=1:n1
 Xi=(1:100)*12*pi+5*rand;Xi=Xi(Xi<len);
 meanISI=mean(diff(Xi));
 Xi = LowPassSmoother(Xi,len,dt,meanISI/2);
 X=[X;Xi];
end
%%% PCA Analysis %%%

%%% Delays + PCA %%%
X2 = [X(:,9:end); X(:,5:end-4); X(:,1:end-8)];
X2=dec(X2,300);
[Xpca,s] = PCA(X2,4);
figure(12);
plot3(Xpca(1,:),Xpca(2,:),Xpca(3,:)+Xpca(4,:),'.','MarkerSize',20);pause
title('PCA (after 3 delays)','Fontsize',18);
Q=[1 1 0 0;1 -1 0 0;0 0 1 1;0 0 1 -1];Q=Q/sqrt(2);
%Xpca(3:6,:)=Q*Xpca(3:6,:);
%[L,diam]=MPlex(Xpca([1 2 3 4 5 6],:),8);
[L,diam]=MPlex(Xpca,8);
