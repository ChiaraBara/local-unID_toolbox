%%   Computes Local Information Dynamics analytically for a stationary mvar(p) process:
%
%   INPUT: 
%   time series Y (one time series in column)
%   ret structure with
%   S_Yn  -   variance of the process

%   OUTPUT:
%   global entropy Hy and local entropy hy

function [Hy, hy] = local_unID_Hlin(Y,ret)

N=size(Y,1);
S_Yn=ret.S_Yn;

% global measures
Hy=0.5*log(2*pi*S_Yn);

% local measures
for n=q+1:N
    yn=Y(n);    % present
    hy(n)=0.5*(yn^2)/S_Yn;
      
end

hy = hy+Hy;

end
