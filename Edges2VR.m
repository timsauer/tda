function [Cd,diam,D]=Edges2VR(numpts,edges,edgediam,tophom)
%numpts is number of 0-simplices
%edges is 2xe matrix of e 1-simplices
S=numpts;Big=160000;tophom2=tophom+2;
C=zeros(tophom2,Big);D=C; % column i of C are indices of points composing simplex i
C(1,1:numpts)=1:numpts;  % column i of D are indices of faces of simplex i
E=zeros(30,Big); % column i of E are indices of simplices for which i is a face
Cd=zeros(1,Big); % Cd(i) is dim of simplex i
Ed=zeros(1,Big); % Ed(i) is length of nonzeros in column i of E 
diam=zeros(1,Big); % diam(i) is diameter of simplex i
[~,numedges]=size(edges);
S=numpts+1;
for i=1:numedges  % while more edges to add
  if 100*round(i/100)==i
    disp([num2str(i) ' of ' num2str(numedges) ' edges'])
  end
  C(1:2,S)=sort(edges(:,i));
  D(1:2,S)=sort(edges(:,i));diam(S)=edgediam(i);
  Cd(S)=1;
  Ed(C(1,S))=Ed(C(1,S))+1;E(Ed(C(1,S)),C(1,S))=S;
  Ed(C(2,S))=Ed(C(2,S))+1;E(Ed(C(2,S)),C(2,S))=S;
  T=S;  %S = number of processed simplices; T = number of simplices found
  while S <= T   % while more unprocessed simplices exist
   if (Cd(S) <= tophom)
    p=Cd(S)+1;  % number of points in simplex
    newsimplex=C(1:p,S); % grab next simplex
    %search for p-simplexes that newsimplex is the final face of
    match=zeros(p,1);cand=[];
    for j=1:p
      fj=D(j,S); %fj=index of jth face
      face=C(1:p-1,fj);% need a common pt that forms p-1 simplex with all p faces
      numk(j)=Ed(fj);cand=union(cand,E(1:numk(j),fj));
      g=find(E(1:numk(j),fj)<=S);numg(j)=length(g); %ignore containing simplices newer than S
      for k=1:numg(j)                               %will get them later
        match(j,k)=setxor(face,C(1:p,E(g(k),fj)));
      end
    end
    inter=match(1,1:numg(1));
    for j=2:p
      inter=intersect(inter,match(j,1:numg(j)));
    end % Find any numbers occuring in all rows of match
     for kk=1:length(inter) % if length(inter)>0 a new simplex has been found
      newsimplex1=sort([inter(kk);newsimplex]); 
      C(1:p+1,T+1)=newsimplex1;
      T=T+1;
      for j=1:p+1
        face=[newsimplex1(1:j-1); newsimplex1(j+1:p+1)];
        [~,~,f]=intersect(face',C(1:p,cand)','rows');
        D(j,T)=cand(f);
        Ed(cand(f))=Ed(cand(f))+1;
        E(Ed(cand(f)),cand(f))=T;
      end
      D(1:p+1,T)=sort(D(1:p+1,T));diam(T)=max(diam(D(1:p+1,T)));
      Cd(T)=p;
     end
   end
   S=S+1;
  end
end
D=D(:,1:T);
Cd=Cd(1:T);
diam=diam(1:T);

      
      
      

      
      
    
    