clear; clc;% close all

addpath([pwd,'\functions\'])
%% Analysis of Time-specific information storage during sleep apnea

load('apnea.mat'); % x=RR, y=RESP
RESP=zscore(y,0,1);

% definition of apnea and respiratory phases
ats=[43 121 202 289 380 488 547 618 681 778 875 972 1040 1119]; %apnea tstart
ate=[82 168 262 355 436 521 579 655 752 841 928 1006 1089 1172];  %apnea tend
ats1=[1  43 83  121 169 202 263 289 356 380 437 488 522 547 580 618 656 681 753 778 842 875 929 972  1007   1040 1090 1119 1173]; 
ate1=[42 82 120 168 201 262 288 355 379 436 487 521 546 579 617 655 680 752 777 841 874 928 971 1006 1039 1089 1118 1172 1200];  

% local IS - linear estimation
q=10; pmax = 20;
[pottaic,pottmdl,aic,mdl] = mos_idMVAR(RESP',pmax,0);
[eAm_R,eSu_R]=idMVAR(RESP',pottmdl,0);
lvarx = local_unID_VAR(eAm_R,eSu_R,q);
out_lin = local_unID_lin(RESP,lvarx);
is_lin = out_lin.sy_Y';

% local IS - knn estimation
k=10; m = 3;
V = [ones(m,1), [1:m]'];
out_knn=local_unID_knn(zscore(RESP),V,k);
is_knn = out_knn.sy_Y;
is_knn = [nan(1,m), is_knn];

%% Figure

lw=1.1; yminy=-4; ymaxy=6; yminx=-2.5; ymaxx=3;
col_as=[0.5 0.5 0.5];
col_ae=[0.8 0.2 0.2];
col1 = [239 111 108]/255;
col2 = [70 87 117]/255;
colorsh=[0.88 0.88 0.88];

t=(1:1:length(RESP))';
figure;
subplot(3,1,1);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[yminy ymaxy ymaxy yminy]',colorsh); hold on
    set(hu,'edgecolor','white');
end
plot(y,'color','k','linewidth',lw); ylim([yminy ymaxy]);
ylabel('RESP')
set(gca, 'YTick', [])

% local IS
sk=2*q;
t=(sk:1:length(y));
ymin=-11; ymax=10;
subplot(3,1,2);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[ymin ymax ymax ymin]',colorsh); hold on
    set(hu,'edgecolor','white');
end
plot(1:1200,is_lin,'color',col1,'linewidth',lw);
ylim([ymin ymax])
ylabel('IS_{lin} [nats]')

% IS knn
ymin=-3; ymax=3;
subplot(3,1,3);
for i=1:length(ats)
    hu=fill([ats(i) ats(i) ate(i) ate(i)]',[ymin ymax ymax ymin]',colorsh); hold on
    set(hu,'edgecolor','white');
end
plot(1:1200,is_knn,'color',col2,'linewidth',lw);
ylim([ymin ymax])
ylabel(' IS_{knn} [nats]')

