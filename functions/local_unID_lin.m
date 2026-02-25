%%   Computes Local Information Dynamics analytically for a stationary mvar(p) process:
%
%   INPUT: 
%   time series Y (one time series in column)
%   ret structure with
%   S_Yn  -   variance of the process
%   S_Ynq  -   covariance matrix of the past of the process 
%   S_YnYnq  -   covariance matrix of joint present and past of the process
%
%   OUTPUT:
%   out structure with local Information Storage

function out = local_unID_lin(Y,ret)

S_Yn=ret.S_Yn;
S_Ynq=ret.S_Ynq;
S_YnYnq=ret.S_YnYnq;

q=size(S_Ynq,1);
N=size(Y,1);

iS_Ynq=inv(S_Ynq);
iS_YnYnq=inv(S_YnYnq);

% global measures
Hy=0.5*log(2*pi*S_Yn);
HY=0.5*log((2*pi*eps(1))^q*det(S_Ynq));
HyY=0.5*log((2*pi*eps(1))^(q+1)*det(S_YnYnq));
Hy_Y=0.5*log((2*pi*eps(1))*det(S_YnYnq)/det(S_Ynq));
%Hy_Y=HyY-HY;
Sy_Y=0.5*log((S_Yn*det(S_Ynq))/det(S_YnYnq));

% local measures
for n=q+1:N
    yn=Y(n);    % present
    ynq=Y(n-1:-1:n-q)'; % past
    
    hy(n)=0.5*(yn^2)/S_Yn;
    hY(n)=0.5*ynq*iS_Ynq*ynq';
    hyY(n)=0.5*[yn ynq]*iS_YnYnq*[yn ynq]';
    hy_Y(n)=0.5*[yn ynq]*iS_YnYnq*[yn ynq]'-0.5*ynq*iS_Ynq*ynq';
%     hy_Y=hyY-hY;
    sy_Y(n)=0.5*(yn^2)/S_Yn + 0.5*ynq*iS_Ynq*ynq'-0.5*[yn ynq]*iS_YnYnq*[yn ynq]';     
end
hy = hy+Hy;
hY = hY+HY;
hyY = hyY+HyY;
hy_Y = hy_Y+Hy_Y;
sy_Y = sy_Y+Sy_Y;

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
