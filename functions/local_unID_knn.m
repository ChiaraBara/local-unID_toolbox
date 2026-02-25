%%   Computes Local Information Dynamics through the nearest-neighbor estimation approach
%
%   INPUT: 
%   time series Y (one time series in column)
%   V assigned embedding vector
%   k number of neighbors
%
%   OUTPUT:
%   out structure with local Information Storage

function out=local_unID_knn(Y,V,k,metric)

if ~exist('metric','var'), metric='chebychev'; end

%% form the observation matrices
B=local_unID_buildvectors(Y,V);
A=B(:,2:end);
N=size(B,1);
m=size(A,2);

tmp=V(:,1);

i_Y= tmp==1;

M_y=B(:,1);
M_Y=A(:,i_Y);
M_yY=[M_y M_Y];

%% kNN analysis
%%% neighbor search in space of higher dimension
[~, distances] = knnsearch(M_yY,M_yY,'K',k+1,'Distance',metric);
dd = distances(:,end);

%%% range searches in subspaces of lower dimension - M_y
[~, distance_y] =  knnsearch(M_y,M_y,'K',N,'Distance',metric);
count_y = distance_y(:,2:end) < dd;
count_y = max(k-1, sum(count_y,2));

%%% range searches in subspaces of lower dimension - M_Y
[~, distance_Y] =  knnsearch(M_Y,M_Y,'K',N,'Distance',metric);
count_Y = distance_Y(:,2:end) < dd;
count_Y = max(k-1, sum(count_Y,2));


%% compute global and local information dynamic measures

% information storage
dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0
Hy= psi(N)-(1/N)*( sum(psi(count_y+1)))+mean(log(dd2));
HY= psi(N)-(1/N)*( sum(psi(count_Y+1)))+m*mean(log(dd2));
HyY= psi(N)-psi(k)+(m+1)*mean(log(dd2));
Hy_Y= - psi(k) + (1/N)*( sum(psi(count_Y+1))) + mean(log(dd2)) ;
Sy_Y = psi(N) + psi(k) -(1/N)*( sum(psi(count_y+1)) + sum(psi(count_Y+1)));
for i = 1:N
    hy(i) = psi(N) - psi(count_y(i)+1) + log(dd2(i));
    hY(i) = psi(N) - psi(count_Y(i)+1) + m*log(dd2(i));
    hyY(i) = psi(N) - psi(k) + (m+1)*log(dd2(i));
    hy_Y(i) = - psi(k) + psi(count_Y(i)+1) + log(dd2(i));
    sy_Y(i) = psi(N) + psi(k) - psi(count_y(i)+1) - psi(count_Y(i)+1);
end

%%% OUTPUT
out.Hy_Y=Hy_Y;
out.Hy=Hy;
out.HY=HY;
out.HyY=HyY;
out.Sy_Y=Sy_Y;

out.hy_Y=hy_Y;
out.hy=hy;
out.hY=hY;
out.hyY=hyY;
out.sy_Y=sy_Y;

end









