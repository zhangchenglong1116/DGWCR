
%¡°Medical Hyperspectral Band Selection Using Ranking Based on Data Gravitation and Weak Correlation¡±

% Input:
% pre     - Raw image data(Preferably normalized), size=A*B, note: A is the number of pixels, B is the number of bands.
% bandnum - Number of band requests
% k       - Powers of the similarity matrix

% Output:
% C       - Band index

function [C] =DGWCR(pre,bandnum,k)
[pix_num,H]=size(pre);                             % H is the number of bands.

%% Information entropy                
x=1;ratio=x/100;                                   % Calculate the information entropy by selecting x % of the pixels evenly.
randNumber = floor(pix_num*ratio);  
idx10=ceil((1:randNumber)*(1/ratio));
pree=pre(idx10,:);
for a =1:H
 bd=pree(:,a);
 rk=unique(bd);
 rn=histc(bd,rk);
 pa(a)=sum(-(rn/pix_num).*log(rn/pix_num));        % Information entropy of bands from Eq.6.
end                   
pa=(pa-min(pa))/(max(pa)-min(pa));                 % Normalized information entropy by Eq.7.

%% Similarity matrix
D=Get_Distance(pre);                               % Distance matrix,Eq.5
D=D/sqrt(pix_num);                                 % Normalized distance matrix
DJ=tril(D);DJ(DJ==0)=[];FOS=sort(DJ); IDX=round((0.02)*length(FOS));c=FOS(IDX);    % The determination of the half-wave peak c of the Gaussian function, with reference to the literature 
                                                                                   % Reference "A Novel Ranking-Based Clustering Approach for Hyperspectral Band Selection" Eq.12
S=exp(-(D.^2)/(c^2));                              % Eq.4.Reference "A Novel Ranking-Based Clustering Approach for Hyperspectral Band Selection" Eq.9

%% CCE search cluster center
SK=S^k;
[~,v]=max(SK);
G=1:H;
CZ=find(v==G);                                     % Eq.1.Reference "Clustering by connection center evolution" Eq.1

%% If there is only one cluster center then sort by self-connectivity
if length(CZ)==1
  [~,Rank]=sort(diag(SK),'descend');
  C=Rank(1:bandnum);
else
    
%% Each point (band) is divided into clusters
i=1;it=1;                                          
while(i~=it+1) 
dense=diag(SK);                                    % self-connectivity
DIAG=dense(CZ);[DENSE,VV]=sort(DIAG,'descend');
CZ=CZ(VV);
cd=SK(CZ,:); Pp=cd./DENSE;
[~,LAB]=max(Pp);class=1:length(CZ);
[class_NUM]=histc(LAB,class);
label=CZ(LAB);                                     
nc=CZ(class_NUM==1);                               % Isolated Point
Fn=CZ(class_NUM>1);                                % Non-isolated point clustering center
CZ=Fn;
i=i+1;
end
mr=max(class_NUM);
RANK=zeros(length(Fn),mr);
des=zeros(length(Fn),1);

%% Intercluster Ranking and Selection of Most Representative Bands based on Data Gravitation
for m =1:length(Fn)
    class_idx=find(label==Fn(m));
    class_idx(class_idx==Fn(m))=[];
    class_idx=[Fn(m) class_idx];
    DEN=pa(class_idx);  
    LD=D(class_idx,class_idx);                     % Local distance matrix from Eq.12
    DistS_1=(LD.*LD).^(-1);
    DistS_1(LD==0)=0;
    Force=sum(DEN.*DistS_1.*DEN');                 % Composite force from Eq.8
    Force=Force/length(Force);
    [sde,VN]=sort(Force,'descend');                % Eq.9
    des(m)=sde(1);
    class_idx=class_idx(VN); 
    RANK(m,1:length(class_idx))=class_idx;         
end

%% Entropy-containing Similarity Matrix.
IF=exp(1./pa);                                     % Eq.11
IS=IF.*S.*IF';                                     % Eq.10

%% Ranking of Other Bands Based on Weak Correlation
for l=1:size(RANK,1)
    HR=RANK(l,:);
    HR(HR==0)=[];
    cs=IS(HR,HR);
    f_nu=1;
    lf=1;
    NR=1:length(HR);
    while(lf<(length(HR)-1))
      sy=setdiff(NR,f_nu); 
      ss=ones(1, length(sy))*10000;
      for t=1:length(sy)
         TM=sy(t);
         js=cs([f_nu TM],[f_nu TM]);
         ss(t)=norm(js,'fro');                     % Eq.14
      end    
    [~,vo]=min(ss);
    vm=sy(vo);                                     % Eq.13
    f_nu=[f_nu vm];
    lf=length(f_nu);
    end
    zone=HR([f_nu setdiff(NR,f_nu)]);
    RANK(l,1:length(HR))=zone;
end

%% S-Shaped Ranking and Band Selection
RANK(RANK==0)=[];                                  % Eq.15 
RANK=[RANK nc];                                    
C=RANK(1:bandnum);
end
end




