function x=dec(x,N)
% Nearest neighbor decimation routine
% Delete each point that is a k-nearest neighbor of more than one point,
% and iterate. Here k=1
s=size(x,2);
while s>N
  w=1:s;
  [ds,inds]=knnCPU(x',x',2);
  h=inds(:,2);u=unique(h);v=u(histcounts(h,u)>1);
  if isempty(v) 
    break; 
  end
  w=setdiff(w,v);
% subplot(4,2,i);plot(x(1,:),x(2,:),'k.',x(1,v),x(2,v),'r.','MarkerSize',30);
  x=x(:,w);
  s=size(x,2);
end