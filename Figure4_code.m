%% Code for figure 4

% Bienzobas Montavez et al. (TurbIFA, submitted to Paleoceanography and
% Paleoclimatology)
%
% --> The sediment mixing-model builds upon SEAMUS software 
% --> The noise estimation (q-q space) and q-q graphical output builds on QUANTIFA software
%
% Source SEAMUS:
% https://github.com/bryanlougheed/seamus
% Lougheed, B. C. (2020). Seamus (v1. 20): a δ14 c-enabled, single-specimen sediment
% accumulation simulator. Geoscientific Model Development, 13 (1), 155–168
%
% Source QUANTIFA:
% https://github.com/rh-glaubke/QUANTIFA
% https://doi.org/10.5281/zenodo.7775163
% https://doi.org/10.1029/2020PA004065
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 25-Dec-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data and declare variables
load('Fig2_4_Tseries.mat','Fig4_Tseries')
full_series_fig4=Fig4_Tseries;

l=length(full_series_fig4);
months_kyrs=12*1000;
num_p=60;                    % number of foraminifera picked
anerr=0;                     % analytical uncertainty 
exerr=0;                     % extra gaussian noise from e.g., calibratiom (1σ,ºC)          
num_q=50;                    % number of quantiles for q-q plot
mc=1000;                     % mc picking times (100 times just to increase speed figures 2-4)                                                   % TurBIFA mc=1000;
run=1;                       % runs of sediment mixing model
X=ones(num_p,1);             % toy data just to get the length of X IFA population (Non-bioturbated)
Y=ones(num_p,1);             % toy data just to get the length of Y IFA population (Bioturbated) 
fig_24=1;

%% Sediment Mixing model panels a,b,c sed.rate 10cm/kyr BD=2cm
sed_r=10;                      % sedimentation rate
BD=2;                          % bioturbation depth (cm)
full_series1=flip(full_series_fig4);
bio_series1=repmat(full_series1,1,run);
z=flip(1:1:(length(bio_series1)/months_kyrs)*sed_r); % pseudo-core
n_months=(months_kyrs/sed_r) ; 
nn_months=repelem(z,n_months);
depths_ind=nn_months;
age=l:-1:1;

% Sediment-mixing
for i = 1:1:length(z)
    for r=1:run
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD);
        a = ind(randperm(length(ind))); % random mixing
        bio_series1(ind,r)=bio_series1(a,r); % bioturbated full_time series
    end
end


full_series1_res=reshape(full_series1',n_months,length(z));
bio_series1_res=reshape(bio_series1,n_months,length(z),run);

% pick data

pos_1= 171;                                          % index discrete depth altered time slice (+2ºC)                                                                                        % position first discrete depth time series 1
bio_series1_bp=squeeze(bio_series1_res(:,pos_1,:));  % bioturbated time series 
full_series1_bp=full_series1_res(:,pos_1);           % non-bioturbated time series

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series1_bp,bio_series1_bp,run,fig_24);

serrX1=pick_output.Xqq_err;
serrY1=pick_output.Yqq_err;
mc_x1=pick_output.mc_x;
mc_y1=pick_output.mc_y;
qx1=pick_output.qqx;
qy1=pick_output.qqy;

% calculate kernels
p1b=zeros(length(mc_x1), 100);
v1b=zeros(length(mc_x1),100);

p1a=zeros(length(mc_y1), 100);
v1a=zeros(length(mc_y1),100);

% before bio
for i=1:1:length(mc_x1)   
   [p1b(i,:),v1b(i,:)]= ksdensity(mc_x1(:,i));
end

% after bio
for i=1:length(mc_y1)
     [p1a(i,:),v1a(i,:)]= ksdensity(mc_y1(:,i));
end



%%  Sediment Mixing model panels d,e,f sed.rate 2cm/kyr BD=10cm
full_series2=full_series1;
bio_series2=repmat(full_series2,1,run);
sed_r=2;
BD=10;
z=flip(1:1:(length(bio_series2)/months_kyrs)*sed_r);
n_months=(months_kyrs/sed_r); 
nn_months=repelem(z,n_months);
depths_ind=nn_months;

% Sediment-mixing
for i = 1:1:length(z)
    for r=1:run
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD);
        a = ind(randperm(length(ind))); % random mixing
        bio_series2(ind,r)=bio_series2(a,r); % bioturbated full_time series
    end
end

full_series2_res=reshape(full_series2',n_months,length(z));
bio_series2_res=reshape(bio_series2,n_months,length(z),run);

% pick

pos_2=35 ;                                           % index discrete depth altered time slice (+2ºC)                                                                                        % position first discrete depth time series 1
bio_series2_bp=squeeze(bio_series2_res(:,pos_2,:));  % bioturbated time series 
full_series2_bp=full_series2_res(:,pos_2);           % non-bioturbated time series

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series2_bp,bio_series2_bp,run,fig_24);

serrX2=pick_output.Xqq_err;
serrY2=pick_output.Yqq_err;
mc_x2=pick_output.mc_x;
mc_y2=pick_output.mc_y;
qx2=pick_output.qqx;
qy2=pick_output.qqy;

% calculate kernels
p2b=zeros(length(mc_x2), 100);
v2b=zeros(length(mc_x2),100);

p2a=zeros(length(mc_y2), 100);
v2a=zeros(length(mc_y2),100);

% before bio
for i=1:1:length(mc_x2)   
   [p2b(i,:),v2b(i,:)]= ksdensity(mc_x2(:,i));
end

% after bio
for i=1:length(mc_y2)
     [p2a(i,:),v2a(i,:)]= ksdensity(mc_y2(:,i));
end


%% Sediment Mixing model panels g h i sed.R 10cm/kyr BD=10cm

full_series3=full_series1;
bio_series3=repmat(full_series3,1,run);
sed_r=10;
BD=10;
z=flip(1:1:(length(bio_series3)/months_kyrs)*sed_r);
n_months=(months_kyrs/sed_r);
nn_months=repelem(z,n_months);
depths_ind=nn_months;

% Sediment-mixing
for i = 1:1:length(z)
    for r=1:run
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+ BD);
        a = ind(randperm(length(ind))); 
        bio_series3(ind,r)=bio_series3(a,r); % bioturbated full_time series
    end
end

full_series3_res=reshape(full_series3',n_months,length(z));
bio_series3_res=reshape(bio_series3,n_months,length(z),run);

% pick

pos_3=171;                                           % index discrete depth altered time slice (+2ºC)                                                                                        % position first discrete depth time series 1
bio_series3_bp=squeeze(bio_series3_res(:,pos_3,:));  % bioturbated time series 
full_series3_bp=full_series3_res(:,pos_3);           % non-bioturbated time series

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series3_bp,bio_series3_bp,run,fig_24);

serrX3=pick_output.Xqq_err;
serrY3=pick_output.Yqq_err;
mc_x3=pick_output.mc_x;
mc_y3=pick_output.mc_y;
qx3=pick_output.qqx;
qy3=pick_output.qqy;

% calculate kernels
p3b=zeros(length(mc_x3),100);
v3b=zeros(length(mc_x3),100);

p3a=zeros(length(mc_y3),100);
v3a=zeros(length(mc_y3),100);

% before bio
for i=1:1:length(mc_x3)   
   [p3b(i,:),v3b(i,:)]= ksdensity(mc_x3(:,i));
end

% after bio
for i=1:length(mc_y3)
     [p3a(i,:),v3a(i,:)]= ksdensity(mc_y3(:,i));
end


%% Sediment Mixing model panels j k l sed.R 2cm/kyr BD=2cm
full_series4=full_series1;
bio_series4=repmat(full_series4,1,run);
sed_r=2;
BD=2;
z=flip(1:1:(length(bio_series4)/months_kyrs)*sed_r);
n_months=(months_kyrs/sed_r) ;  
nn_months=repelem(z,n_months);
depths_ind=nn_months;

% Sediment-mixing
for i = 1:1:length(z)
    for r=1:run
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD);
        a = ind(randperm(length(ind))); 
        bio_series4(ind,r)=bio_series4(a,r); %bioturbated full_time series
    end
end

full_series4_res=reshape(full_series4',n_months,length(z)); % nº months per cm
bio_series4_res=reshape(bio_series4,n_months,length(z),run);

% pick

pos_4=35;                                            % index discrete depth altered time slice (+2ºC)                                                                                        % position first discrete depth time series 1
bio_series4_bp=squeeze(bio_series4_res(:,pos_4,:));  % bioturbated time series 
full_series4_bp=full_series4_res(:,pos_4);           % non-bioturbated time series

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series4_bp,bio_series4_bp,run,fig_24);

serrX4=pick_output.Xqq_err;
serrY4=pick_output.Yqq_err;
mc_x4=pick_output.mc_x;
mc_y4=pick_output.mc_y;
qx4=pick_output.qqx;
qy4=pick_output.qqy;

% calculate kernels
p4b=zeros(length(mc_x4),100);
v4b=zeros(length(mc_x4),100);

p4a=zeros(length(mc_y4),100);
v4a=zeros(length(mc_y4),100);

% before bio
for i=1:1:length(mc_x4)   
   [p4b(i,:),v4b(i,:)]= ksdensity(mc_x4(:,i));
end

% after bio
for i=1:length(mc_y4)
     [p4a(i,:),v4a(i,:)]= ksdensity(mc_y4(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data

subplot(4,3,1)
plot(flip(age),flip(bio_series1(:,1)),'Color','0.8 0.8 0.8');
hold on;
plot(flip(age),flip(full_series1),'Color','black');
xlabel('month');
ylabel('T (ºC)');
ylim([16 37])
title('10cm/kyr | 2cm')

subplot(4,3,2)

for i=1:1:length(v1b) % before bio
    plot(v1b(i,:),p1b(i,:),'Color','[0, 0, 1, 0.03]')
    hold on
end

hold on
for i=1:1:length(v1a) % after bio 
    plot(v1a(i,:),p1a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end
hold on
%plot(mean(v1b),mean(p1b),'Color','b','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v1a),mean(p1a),'Color','red','LineWidth',2) % mean kernels after bio
xlabel('T (°C)');
ylabel('Prob.Density');
xlim([15 37])
ylim([0 0.3])


subplot(4,3,3)
errorbar(median(qx1,2),median(qy1,2),serrY1(1,:),serrY1(2,:),'linestyle','none','color','b','capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx1,2),median(qy1,2),serrX1(1,:),serrX1(2,:),'horiz','color','b','linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx1,2),median(qy1,2),'kd-','MarkerSize',5,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (°C)');
ylabel('After Bio (°C)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([18 35])
ylim([18 35])

subplot(4,3,4)
plot(flip(age),flip(bio_series2(:,1)),'Color','0.8 0.8 0.8');
hold on;
plot(flip(age),flip(full_series2),'Color','black');
xlabel('month');
ylabel('T (ºC)');
ylim([16 37])
title('2cm/kyr | 10cm')

subplot(4,3,5)
for i=1:1:length(v2b) % before bio
    plot(v2b(i,:),p2b(i,:),'Color','[0, 0, 1, 0.03]')
    hold on
end

hold on
for i=1:1:length(v2a) % after bio 
    plot(v2a(i,:),p2a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on
%plot(mean(v2b),mean(p2b),'Color','b','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v2a),mean(p2a),'Color','red','LineWidth',2) % mean kernels after bio
xlabel('T (°C)');
ylabel('Prob.Density')
xlim([15 37])
ylim([0 0.4])


subplot(4,3,6)

errorbar(median(qx2,2),median(qy2,2),serrY2(1,:),serrY2(2,:),'linestyle','none','color','b','capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx2,2),median(qy2,2),serrX2(1,:),serrX2(2,:),'horiz','color','b','linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx2,2),median(qy2,2),'kd-','MarkerSize',5,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (°C)');
ylabel('After Bio (°C)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([18 35])
ylim([18 35])

subplot(4,3,7)
plot(flip(age),flip(bio_series3(:,1)),'Color','0.8 0.8 0.8');
hold on;
plot(flip(age),flip(full_series3),'Color','black');
xlabel('month');
ylabel('T (ºC)');
ylim([16 37])
title('10cm/kyr | 10cm')

subplot(4,3,8)

for i=1:1:length(v3b) % before bio
    plot(v3b(i,:),p3b(i,:),'Color','[0, 0, 1, 0.03]')
    hold on
end

hold on
for i=1:1:length(v3a) % after bio 
    plot(v3a(i,:),p3a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on
%plot(mean(v3b),mean(p3b),'Color','b','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v3a),mean(p3a),'Color','red','LineWidth',2) % mean kernels after bio
xlabel('T (ºC)');
ylabel('Prob.Density')
xlim([15 37])
ylim([0 0.3])

subplot(4,3,9)
errorbar(median(qx3,2),median(qy3,2),serrY3(1,:),serrY3(2,:),'linestyle','none','color','b','capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx3,2),median(qy3,2),serrX3(1,:),serrX3(2,:),'horiz','color','b','linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx3,2),median(qy3,2),'kd-','MarkerSize',5,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (ºC)');
ylabel('After Bio (ºC)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([18 35])
ylim([18 35])

subplot(4,3,10)
plot(flip(age),flip(bio_series4(:,1)),'Color','0.8 0.8 0.8');
hold on;
plot(flip(age),flip(full_series4),'Color','black');
xlabel('month');
ylabel('T (ºC)');
ylim([16 37])
title('2cm/kyr | 2cm')

subplot(4,3,11)

for i=1:1:length(v4b) % before bio
    plot(v4b(i,:),p4b(i,:),'Color','[0, 0, 1, 0.03]')
    hold on
end

hold on
for i=1:1:length(v4b) % after bio 
    plot(v4a(i,:),p4a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on
%plot(mean(v4b),mean(p4b),'Color','b','LineWidth',2)   % mean kernels before bio
hold on
%plot(mean(v4a),mean(p4a),'Color','red','LineWidth',2) % mean kernels after bio
xlabel('T (ºC)');
ylabel('Prob.Density')
xlim([15 37])
ylim([0 0.3])

subplot(4,3,12)

errorbar(median(qx4,2),median(qy4,2),serrY4(1,:),serrY4(2,:),'linestyle','none','color','b','capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx4,2),median(qy4,2),serrX4(1,:),serrX4(2,:),'horiz','color','b','linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx4,2),median(qy4,2),'kd-','MarkerSize',5,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (ºC)');
ylabel('After Bio (ºC)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([18 35])
ylim([18 35])


annotation('textbox', [0.08,0.84,0.21,0.15], 'String', 'a)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.84,0.21,0.15], 'String', 'b)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.84,0.21,0.15], 'String', 'c)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.62,0.21,0.15], 'String', 'd)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.62,0.21,0.15], 'String', 'e)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.62,0.21,0.15], 'String', 'f)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.40,0.21,0.15], 'String', 'g)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.40,0.21,0.15], 'String', 'h)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.40,0.21,0.15], 'String', 'i)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.18,0.21,0.15], 'String', 'j)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.18,0.21,0.15], 'String', 'k)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.18,0.21,0.15], 'String', 'l)', 'FontSize', 9, 'EdgeColor', 'none');

clear;