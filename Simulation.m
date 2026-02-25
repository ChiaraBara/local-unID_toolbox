%% Test Simulation
% Time-specific information storage with innovation variance and AR coefficient 
% which varying as a periodic square waveform 
close all;clear all;clc

addpath([pwd,'\functions\']);

%% generation of time variant innovation variance and AR coefficient

rng(1)
cmax      =  0.3;    % max amplitude
fs        =  1000;    % sampling frequency
f_osc     =  0.05;    % frequency of oscillation 
p         =    1;    % model order        

%% Construction of non-stationary VAR time series

nobs = 60000;                        % 5 s of observation (data samples)
k = 1:nobs;                         % vector of time steps
t = (k)/fs;                         % vector of times
DC=50;                              % duty-cycle
c = -cmax*square(2*pi*f_osc*t,DC);  % causal coefficient varying as square periodic waveform
Su_T(1,1,:)=c+cmax*2;               % theoretical vector of innovation variance 
A_T(1,1,:) = c+cmax*2;              % theoretical AR coefficient               
% Data generation
rng(1) % always get the same seed to fix the realization
Y = var_nonstat(A_T,Su_T,1);
% theoretical IS value
ret_t=tv_IS(A_T,Su_T);

%% local IS - linear
Y=squeeze(Y)';
[eAm,eSu]=idMVAR(Y,p,0);
lvar = local_unID_VAR(eAm,eSu,10);
out_lin = local_unID_lin(Y',lvar);
is_lin = out_lin.sy_Y;

LIM=[-15, 15];
figure;
xlim([0.2 60])
ylim(LIM)
D=abs(diff(c));
[ind]=find(D==max(D));
IND=zeros(1,length(ind)+1);
IND(1)=1;
IND(2:end)=ind;
MEAN=[12];
VAR=[-12];
pts=200;
for oo=1:length(IND)-1
    M=nanmean(is_lin(IND(oo):IND(oo+1)));
    V=nanvar(is_lin(IND(oo):IND(oo+1)));
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),MEAN,sprintf('%4.3f',round(M,3,'Significant')) ,'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),VAR,sprintf('(%4.3f)',round(V,3,'Significant')),'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
end
hold on
for pp=1:numel(ind)
    xline(t(ind(pp)),'--k')
end
hold on
plot(t(pts:end),is_lin(pts:end),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
hold on
plot(t(pts:end),[ret_t.IS(pts:end,1)],'k','LineWidth',1);
ylabel('IS_{lin} [nats]')
xlabel('Time [s]')
set(gcf,'units','centimeters','position',[0,0,12,12])

%% local IS - KNN
k = 10;
m = 1;

Y=squeeze(Y)';
V = [ones(m,1), [1:m]'];
out_knn=local_unID_knn_mex(Y,V,k);
is_knn = out_knn.sy_Y;
is_knn = [nan(1,m), is_knn];

LIM=[-15, 15];
figure;
xlim([0.2 60])
ylim(LIM)
D=abs(diff(c));
[ind]=find(D==max(D));
IND=zeros(1,length(ind)+1);
IND(1)=1;
IND(2:end)=ind;
MEAN=[12];
VAR=[-12];
pts=200;
for oo=1:length(IND)-1
    M=nanmean(is_knn(IND(oo):IND(oo+1)));
    V=nanvar(is_knn(IND(oo):IND(oo+1)));
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),MEAN,sprintf('%4.3f',round(M,3,'Significant')) ,'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
    text(t(round(IND(oo+1)-((IND(oo+1)-IND(oo))/2))),VAR,sprintf('(%4.3f)',round(V,3,'Significant')),'FontSize',...
        10,'Color','k','FontName','TimesNewRoman','HorizontalAlignment','center');
end
hold on
for pp=1:numel(ind)
    xline(t(ind(pp)),'--k')
end
hold on
plot(t(pts:end),is_knn(pts:end),'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
hold on
plot(t(pts:end),[ret_t.IS(pts:end,1)],'k','LineWidth',1);
ylabel('IS_{knn} [nats]')
xlabel('Time [s]')
set(gcf,'units','centimeters','position',[0,0,12,12])


figure; 
plot(t,squeeze(A_T),'b','LineWidth',1);
hold on; plot(t,squeeze(Su_T),'--r','LineWidth',1);
for pp=1:numel(ind)
    xline(t(ind(pp)),'--k')
end
legend({'a_n','\sigma^2_n'})
ylim([0 1]);
