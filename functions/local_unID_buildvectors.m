%% form observation matrix 
%%% INPUT
% Y: input multiple time series, dimension N*M
% V: list of candidates, dimension Nc*2, Nc is number of candidates; 1st column: index of the signal; 2nd column: index of the lag
%%% OUTPUT
% B: complete matrix with added the current samples as first column

function B=local_unID_buildvectors(Y,V)

if isempty(V) %if no conditioning, give back the jth series of Y
    B=Y;
else
    [N,M]=size(Y);
    Nc=size(V,1); % number of candidates

    Lmax=max(V(:,2)); %maximum lag (across all signals)

    A=NaN*ones(N-Lmax,Nc);
    for n=Lmax+1:N
        for i=1:Nc %fill the i-th row of A
            A(n-Lmax,i)=Y(n-V(i,2),V(i,1)); 
        end
    end
    B=[Y(Lmax+1:N) A]; % add current value
end
