function [L,diam]=MPlex(data,meanEdgesPerNode,k,k2,tophom,delta);
% assume the points are in dim x number of pts variable called data
% k2=number of nbrs in CkNN, must be <= k=number of nbrs found in knnCPU
% tophom=top dimension of homology desired
finddelta=0; %if finddelta = 0, use given delta
if nargin<6
  finddelta=1; %if finddelta = 1, compute a delta below
  if nargin<5
    tophom=3;
    if nargin<4
      k2=5;
      if nargin<3
        k=5;
        if nargin<2
          meanEdgesPerNode=6;
        end
      end
    end
  end
end
close all
[datadim,numpts]=size(data);
[d,inds]=knnCPU(data',data',k);
d=sqrt(d(:,k2));d=diag(1./d);
dist=1000000*ones(numpts,numpts);
for i=1:numpts
  for j=i+1:numpts
    dist(i,j)=norm(data(:,i)-data(:,j),2);
  end
end
dist=d*dist*d;   %CkNN (omit this line for conventional e-ball)
[edgediam,ra]=sort(dist(:));
[I,J]=ind2sub(size(dist),ra);di=round(meanEdgesPerNode*numpts);
edges=[I(1:di) J(1:di)]';edgediam=edgediam(1:di);edgediam=(1:di);
tic;[dimen,diam,bound]=Edges2VR(numpts,edges,edgediam,tophom);toc
[q,p]=size(bound);
time=[ones(1,numpts) (1:(p-numpts))];time=diam;
tic;L=ComputeIntervals(dimen,bound,p,time);toc
maxtime=max(time);
% Betti=PlotBarCodes(L,maxtime+.5,tophom);
% disp(['Betti numbers in increasing order of dim: ']);Betti
while di>0
 Betti=PlotBarCodes(L,di,numpts,tophom);
 disp(['Betti numbers in increasing order of dim: ']);Betti
 figure(9);
if datadim>=3
  plot3(data(1,:),data(2,:),data(3,:),'b.','MarkerSize',25);hold on
  for i=1:di
    plot3([data(1,I(i)) data(1,J(i))],[data(2,I(i)) data(2,J(i))],[data(3,I(i)) data(3,J(i))],'r','LineWidth',2)
  end
  plot3(data(1,:),data(2,:),data(3,:),'b.','MarkerSize',25);hold off
elseif datadim==2
  plot(data(1,:),data(2,:),'b.','MarkerSize',25);hold on
  for i=1:di
    plot([data(1,I(i)) data(1,J(i))],[data(2,I(i)) data(2,J(i))],'r','LineWidth',2)
  end
  plot(data(1,:),data(2,:),'b.','MarkerSize',25);hold off
end
di=input('Show graph with this reduced number of mean edges (stop = 0):');
di=round(numpts*di);
end

