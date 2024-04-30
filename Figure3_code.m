%% Code for figure 3

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
% Author: Natalia Bienzobas Montávez 
%
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 25-Dec-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data and declare variables

load('Fig2_4_Tseries.mat','Fig3_Tseries1','Fig3_Tseries2')
full_series1=Fig3_Tseries1; % 1-kyr time slice mean temperature altered +2deg.
full_series2=Fig3_Tseries2; % 1-kyr time slice mean temperature altered -2deg.

% Declare Data

l=length(full_series1);
num_p=60;                    % number of foraminifera picked
months_kyrs=12*1000;         % number of months per kyr
sed_r=2;                     % sedimentation rate (cm/kyr)
BD=5;                        % SMLD/bioturbation depth (cm)
anerr=0;                     % analytical uncertainty 
exerr = 0;                   % extra gaussian noise from e.g., calibratiom (1σ,ºC) 
num_q=50;                    % number of quantile for q-q plot
mc=1000;                     % mc picking times (100 times just to increase speed figures 2-4)                                                   % TurBIFA mc=1000;
run=1;                       % runs of sediment mixing model
X=ones(num_p,1);             % toy data just to get the length of X IFA population (Non-bioturbated)
Y=ones(num_p,1);             % toy data just to get the length of Y IFA population (Bioturbated) 
fig_24=1;
% Note that in TurBIFA T and d18O models, X and Y variables corresponds to
% real IFA populations

z=flip(1:1:(length(full_series1)/12000)*sed_r);  % pseudo sediment core
n_months=(months_kyrs/sed_r) ;  % nº month per cm
nn_months=repelem(z,n_months);
depths_ind=nn_months;
age=l:-1:1;
age_res=reshape(age,n_months,length(z));  
bio_series1=repmat(full_series1,1,run);
bio_series2=repmat(full_series2,1,run);

%%
for r=1:run
    for i = 1:1:length(z)
         ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD);
         a = ind(randperm(length(ind))); % random mixing
         bio_series1(ind,r)=bio_series1(a,r); %bioturbated series_1
         bio_series2(ind,r)=bio_series2(a,r); %bioturbated series_2
     end
end

bio_series1_res=reshape(bio_series1,n_months,length(z),run);
full_series1_res=reshape(full_series1',n_months,length(z));
bio_series2_res=reshape(bio_series2,n_months,length(z),run);
full_series2_res=reshape(full_series2',n_months,length(z));
%% Resampling ("picking")

% 1. Upper panels
pos_1= 35;                                           % index discrete depth altered time slice (+2ºC)                                                                                        % position first discrete depth time series 1
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
for i=1:1:length(mc_y1)
     [p1a(i,:),v1a(i,:)]= ksdensity(mc_y1(:,i));
end

% 2. Lower panels

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
for i=1:1:length(mc_y2)
     [p2a(i,:),v2a(i,:)]= ksdensity(mc_y2(:,i));
end

%% Plot
subplot(2,3,1)
plot(flip(age),flip(bio_series1(:,1)),'Color','0.8 0.8 0.8'); % plot 1 run of sediment mixing model
hold on
plot(flip(age),flip(full_series1),'Color','black');
hold on
x2 = [min(min(age_res(:,35:36))) max(max(age_res(:,35:36))) max(max(age_res(:,35:36))) min(min(age_res(:,35:36)))];
y2 = [min(min(full_series1_res(:,35:36)))    min(min(full_series1_res(:,35:36)))  max(max(full_series1_res(:,35:36)))  max(max(full_series1_res(:,35:36)))];
patch(x2,y2,[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880])
ylim([20 35])
xlabel('month');
ylabel('T (°C)');


subplot(2,3,2)

for i=1:1:length(v1b) % before bio
    plot(v1b(i,:),p1b(i,:),'Color','[0.4660 0.6740 0.1880 0.03]')
    hold on
end

hold on
for i=1:1:length(v1a) 
    plot(v1a(i,:),p1a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end
hold on
%plot(mean(v1b),mean(p1b),'Color','[0.4660 0.6740 0.1880]','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v1a),mean(p1a),'Color','red','LineWidth',2) % mean kernel after bio
xlabel('T (°C)');
ylabel('Prob.Density');
xlim([15 35])
ylim([0 0.4])

subplot(2,3,3)
errorbar(median(qx1,2),median(qy1,2),serrY1(1,:),serrY1(2,:),'linestyle','none','color',[0.4660 0.6740 0.1880],'capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx1,2),median(qy1,2),serrX1(1,:),serrX1(2,:),'horiz','color',[0.4660 0.6740 0.1880],'linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx1,2),median(qy1,2),'kd-','MarkerSize',6,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (°C)');
ylabel('After Bio (°C)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
xlim([21.5 31.5])
ylim([21.5 31.5])


subplot(2,3,4)
plot(flip(age),flip(bio_series2),'Color','0.8 0.8 0.8');
hold on
plot(flip(age),flip(full_series2),'Color','black');
hold on
x2 = [min(min(age_res(:,35:36))) max(max(age_res(:,35:36))) max(max(age_res(:,35:36))) min(min(age_res(:,35:36)))];
y2 = [min(min(full_series2_res(:,35:36)))   min(min(full_series2_res(:,35:36)))   max(max(full_series2_res(:,35:36)))  max(max(full_series2_res(:,35:36)))];
patch(x2,y2,[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250])
ylim([17 30])
xlabel('month');
ylabel('T (ºC)');

subplot(2,3,5)

for i=1:1:length(v2b) % before bio
    plot(v2b(i,:),p2b(i,:),'Color','[0.9290 0.6940 0.1250 0.03]')
    hold on
end

hold on
for i=1:1:length(v2a) 
    plot(v2a(i,:),p2a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on
%plot(mean(v2b),mean(p2b),'Color','[0.9290 0.6940 0.1250]','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v2a),mean(p2a),'Color','red','LineWidth',2) % mean kernels after bio
xlabel('T (°C)');
ylabel('Prob.Density')
xlim([15 35])
ylim([0 0.4])


subplot(2,3,6)
errorbar(median(qx2,2),median(qy2,2),serrY2(1,:),serrY2(2,:),'linestyle','none','color',[0.9290 0.6940 0.1250],'capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx2,2),median(qy2,2),serrX2(1,:),serrX2(2,:),'horiz','color',[0.9290 0.6940 0.1250],'linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx2,2),median(qy2,2),'kd-','MarkerSize',6,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (°C)');
ylabel('After Bio (°C)');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
xlim([19 30])
ylim([19 30])


annotation('textbox', [0.08,0.67,0.21,0.34], 'String', 'a)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.67,0.21,0.34], 'String', 'b)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.67,0.21,0.34], 'String', 'c)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.19,0.21,0.34], 'String', 'd)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.19,0.21,0.34], 'String', 'e)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.19,0.21,0.34], 'String', 'f)', 'FontSize', 9, 'EdgeColor', 'none');

clear;