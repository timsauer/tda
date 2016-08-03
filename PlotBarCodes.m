function Betti=PlotBarCodes(L,t,numpts,tophom)
maxdim=max(L(4,:));maxdim=min(maxdim,tophom);
maxtime=max(1,max(L(2,:)));
for dim=0:maxdim
  f=find((L(4,:)==dim).*(L(2,:)>L(1,:)));num=length(f);L1=L(:,f);
  [~,r]=sort(L1(1,:),'descend');L1=L1(:,r);
  y=linspace(0,1,num+2);y=y(2:end-1);
  figure(dim+1);hold on
  for i=1:num
    plot([L1(1,i) L1(2,i)+.05]/numpts,[y(i) y(i)],'LineWidth',2)
  end
  hold off
  axis([0 maxtime/numpts 0 1]);grid on;title(['Homology in dim ' num2str(dim)],'FontSize',20)
  xlabel('mean edges per node')
end
Betti=zeros(maxdim+1,1);
for dim=0:maxdim
  f=find((L(4,:)==dim).*(L(2,:)>L(1,:)));num=length(f);L1=L(:,f);
  counter=0;
  for i=1:num
    if L1(1,i)<=t && t<=L1(2,i)
      counter=counter+1;
    end
  end
  Betti(dim+1)=counter;
end
