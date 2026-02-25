%%   Computes Local Information Dynamics through the nearest-neighbor estimation approach
%
%   INPUT: 
%   time series Y (one time series in column)
%   k number of neighbors
%
%   OUTPUT:
%   global entropy Hy and local entropy hy

function [Hy, hy]=local_unID_Hknn(Y,k,metric)

if ~exist('metric','var'), metric='chebychev'; end

N=size(Y,1);

%% kNN analysis
%%% neighbor search in space of higher dimension
[~, distances] = knnsearch(Y,Y,'K',k+1,'Distance',metric);
dd = distances(:,end);

% information storage
dd2=2*dd;
dd2(dd2==0)=[]; % do not accept distance=0

Hy= psi(N)-psi(k)+mean(log(dd2));
for i = 1:N
    hy(i) = psi(N) - psi(k) + log(dd2(i));
end

end









