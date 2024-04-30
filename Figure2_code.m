%% Code for figure 2

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
% Last revision: 24-Dec-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data and declare variables

load('Fig2_4_Tseries.mat','Fig2_Tseries')
full_series=Fig2_Tseries;

num_p=60;                  % number of foraminifera picked
months_kyrs=1000*12;       % number of months per kyr
sed_r=10;                  % sedimentation rate cm/kyr
BD=15;                     % SMLD/bioturbation depth (cm)
anerr=0;                   % analytical uncertainty (1σ,ºC)
exerr = 0;                 % extra gaussian noise from e.g., calibratiom (1σ,ºC) 
num_q=50;                  % number of quantile for q-q plot
mc=1000;                   % mc picking times (100 times just to increase speed figures 2-4)                                                   % TurBIFA mc=1000;
run=1;                     % runs of sediment mixing model
X=ones(num_p,1);           % toy data just to get the length of X IFA population (Non-bioturbated)
Y=ones(num_p,1);           % toy data just to get the length of Y IFA population (Bioturbated) 
fig_24=1;

% Note that in TurBIFA T and d18O models, X and Y variables corresponds to
% real IFA populations

z=flip(1:1:(length(full_series)/12000)*sed_r);  % pseudo sediment core
n_months=(months_kyrs/sed_r) ;                  
nn_months=repelem(z,n_months);
depths_ind=nn_months;
l=length(full_series);
age=l:-1:1;
age_res=reshape(age,n_months,length(z));  % months per cm 

%1 --> Figure 2 Panels a,b,c Bioturbated first 6cm

% Declare data
z1=z(1:15);        % first 6cm (oldest)
full_series1=flip(full_series);  
bio_series1=repmat(full_series1,1,run);

% Sediment-mixing model first 6cm
for r= 1:run  % repeat mixing run times 
    for i = 1:1:length(z1)
        ind= find(depths_ind >= z1(i) & depths_ind <z1(i)+BD);
        a = ind(randperm(length(ind))); 
        bio_series1(ind,r)=bio_series1(a,r); 
    end
end


bio_series1_res=reshape(bio_series1,n_months,length(z),run);
full_series1_res=reshape(full_series1',n_months,length(z));

% 2 --> Figure 2 Panels d,e,f Bioturbated until the depth where the
% amplitude of variability altered (+2x) 

z2=z(1:115);  % from core bottom to first 6cm of altered (+2x) 1kyr time slice  
full_series2=flip(full_series); 
bio_series2=repmat(full_series2,1,run);
 
% Sediment-mixing model 

for r= 1:run  % repeat mixing run times 
    for i = 1:1:length(z2)
        ind= find(depths_ind >= z2(i) & depths_ind <z2(i)+BD);
        a = ind(randperm(length(ind))); 
        bio_series2(ind,r)=bio_series2(a,r); 
    end
end   


bio_series2_res=reshape(bio_series2,n_months,length(z),run);
full_series2_res=reshape(full_series2',n_months,length(z));


% 3 --> Figure 3 Panels g,h,i Bioturbated until the depth where the
% amplitud of variability altered (-2x)

z3=z(1:215);  % from core bottom to first 6cm of altered -2x 1kyr time slice     
full_series3=flip(full_series);
bio_series3=repmat(full_series3,1,run);
 
% Sediment-mixing model 

for r= 1:run  % repeat mixing run times 
    for i = 1:1:length(z3)
          ind= find(depths_ind >= z3(i) & depths_ind <z3(i)+BD);
          a = ind(randperm(length(ind))); 
          bio_series3(ind,run)=bio_series3(a,run); 
    end
end
   

bio_series3_res=reshape(bio_series3,n_months,length(z),run);
full_series3_res=reshape(full_series3',n_months,length(z));

%% Resampling ("picking")
% Panels a,b,c

% Lets say that we pick forams over a discrete depth inside the active SMLD (first 6cm)
pos_1= 1;           % index position  discrete depth                                                                                         % position first discrete depth time series 1
bio_series1_bp=squeeze(bio_series1_res(:,pos_1,:));  % bioturbated time series 1 before picking
full_series1_bp=full_series1_res(:,pos_1);            

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series1_bp,bio_series1_bp,run,fig_24);


serrX1=pick_output.Xqq_err;
serrY1=pick_output.Yqq_err;
mc_x1=pick_output.mc_x;
mc_y1=pick_output.mc_y;
qx1=pick_output.qqx;
qy1=pick_output.qqy;

%Panels d,e,f

pos_2= 101;           % position first discrete depth inside +2x time slice
bio_series2_bp=squeeze(bio_series2_res(:,pos_2,:));  % bioturbated time series 1 before picking
full_series2_bp=full_series2_res(:,pos_2);            

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series2_bp,bio_series2_bp,run,fig_24);

serrX2=pick_output.Xqq_err;
serrY2=pick_output.Yqq_err;
mc_x2=pick_output.mc_x;
mc_y2=pick_output.mc_y;
qx2=pick_output.qqx;
qy2=pick_output.qqy;


%Panels g,h,i

pos_3= 201;           % position first discrete depth inside -2x time slice

bio_series3_bp=squeeze(bio_series3_res(:,pos_3,:));  % bioturbated time series 1 before picking
full_series3_bp=full_series3_res(:,pos_3);            

pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,full_series3_bp,bio_series3_bp,run,fig_24);

serrX3=pick_output.Xqq_err;
serrY3=pick_output.Yqq_err;
mc_x3=pick_output.mc_x;
mc_y3=pick_output.mc_y;
qx3=pick_output.qqx;
qy3=pick_output.qqy;


%% Some calculations before final plot

% 1. panel b  Kernel density 

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

% 2. panel d Kenel density
p2b=zeros(length(mc_x2), 100);
v2b=zeros(length(mc_x2),100);

p2a=zeros(length(mc_x2), 100);
v2a=zeros(length(mc_y2),100);

% before bio
for i=1:1:length(mc_x2)   
   [p2b(i,:),v2b(i,:)]= ksdensity(mc_x2(:,i));
end

% after bio
for i=1:length(mc_y2)
   [p2a(i,:),v2a(i,:)]= ksdensity(mc_y2(:,i));
end

% 3. panel d Kenel density
p3b=zeros(length(mc_x3), 100);
v3b=zeros(length(mc_x3),100);

p3a=zeros(length(mc_y3), 100);
v3a=zeros(length(mc_y3),100);

% before bio
for i=1:1:length(mc_x3)   
   [p3b(i,:),v3b(i,:)]= ksdensity(mc_x3(:,i));
end

% after bio
for i=1:length(mc_y3)
   [p3a(i,:),v3a(i,:)]= ksdensity(mc_y3(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot results 


subplot(3,3,1)
plot(flip(age),full_series,'Color','black');
hold on
x2 = [min(min(age_res(:,101:110))) max(max(age_res(:,101:110))) max(max(age_res(:,101:110))) min(min(age_res(:,101:110)))];
y2 = [min(min(full_series2_res(:,101:110)))    min(min(full_series2_res(:,101:110)))  max(max(full_series2_res(:,101:110)))  max(max(full_series2_res(:,101:110)))];
patch(x2,y2,[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880])

hold on
x2 = [min(min(age_res(:,201:210))) max(max(age_res(:,201:210))) max(max(age_res(:,201:210))) min(min(age_res(:,201:210)))];
y2 = [min(min(full_series3_res(:,201:210)))    min(min(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))];
patch(x2,y2,[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250])

hold on
x2 = [min(min(age_res(:,1:15))) max(max(age_res(:,1:15))) max(max(age_res(:,1:15))) min(min(age_res(:,1:15)))];
y2 = [min(min(full_series1_res(:,1:15)))    min(min(full_series1_res(:,1:15)))  max(max(full_series1_res(:,1:15)))  max(max(full_series1_res(:,1:15)))];
patch(x2,y2,'red','EdgeColor','red')

xlabel('Months');
ylabel('Temperature (°C)')
set(gca, 'XAxisLocation', 'top')
ylim([15 35])

subplot(3,3,2)
for i=1:1:length(v1a) % after bio 
    plot(v1a(i,:),p1a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end
hold on

for i=1:1:length(v1b) % before bio
    plot(v1b(i,:),p1b(i,:),'Color','[0, 0, 0.8, 0.03]')
    hold on
end


hold on
%plot(mean(v1b),mean(p1b),'Color','blue','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v1a),mean(p1a),'Color','red','LineWidth',2) % mean kernels after bio

xlim([15 35])
ylim([0 0.45])
xlabel('Temperature (ºC)');
ylabel('Prob.Density')
set(gca, 'XAxisLocation', 'top')


subplot(3,3,3)
errorbar(median(qx1,2),median(qy1,2),serrY1(1,:),serrY1(2,:),'linestyle','none','color','b','capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx1,2),median(qy1,2),serrX1(1,:),serrX1(2,:),'horiz','color','b','linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx1,2),median(qy1,2),'kd-','MarkerSize',6,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (ºC)','Color','b');
ylabel('After Bio (ºC)','Color','r');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
xlim([19 31])
ylim([19 31])
r1.XData=19:31;
r1.YData=19:31;


subplot(3,3,4)
plot(flip(age),full_series,'Color','black');
hold on
data2_bio=(bio_series2(:,1));  
plot(flip(age(1:(n_months*115))),flip(data2_bio(1:(n_months*115))),'Color','0.8 0.8 0.8');
hold on
x2 = [min(min(age_res(:,101:115))) max(max(age_res(:,101:115))) max(max(age_res(:,101:115))) min(min(age_res(:,101:115)))];
y2 = [min(min(full_series2_res(:,101:115)))    min(min(full_series2_res(:,101:115)))  max(max(full_series2_res(:,101:115)))  max(max(full_series2_res(:,101:115)))];
patch(x2,y2,'red','EdgeColor','red')
hold on
x2 = [min(min(age_res(:,101:110))) max(max(age_res(:,101:110))) max(max(age_res(:,101:110))) min(min(age_res(:,101:110)))];
y2 = [min(min(full_series2_res(:,101:110)))    min(min(full_series2_res(:,101:110)))  max(max(full_series2_res(:,101:110)))  max(max(full_series2_res(:,101:110)))];
patch(x2,y2,[0.4660 0.6740 0.1880],'EdgeColor',[0.4660 0.6740 0.1880],'FaceColor','none')
xlabel('Months');
ylabel('Temperature (°C)')
ylim([15 35])

subplot(3,3,5)

for i=1:1:length(v2a) % after bio 
    plot(v2a(i,:),p2a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on

for i=1:1:length(v2b) % before bio
    plot(v2b(i,:),p2b(i,:),'Color','[0.4660 0.6740 0.1880 0.03]')
    hold on
end

hold on
%plot(mean(v2b),mean(p2b),'Color','[0.4660 0.6740 0.1880]','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v2a),mean(p2a),'Color','red','LineWidth',2) % mean kernels after bio

xlim([15 35])
ylim([0 0.27])
xlabel('Temperature (°C)')
ylabel('Prob.Density')


subplot(3,3,6)
errorbar(median(qx2,2),median(qy2,2),serrY2(1,:),serrY2(2,:),'linestyle','none','color',[0.4660 0.6740 0.1880],'capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx2,2),median(qy2,2),serrX2(1,:),serrX2(2,:),'horiz','color',[0.4660 0.6740 0.1880],'linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx2,2),median(qy2,2),'kd-','MarkerSize',6,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (ºC)','Color','[0.4660 0.6740 0.1880]');
ylabel('After Bio (ºC)','Color','r');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([15 35])
ylim([15 35])
r1.XData=15:35;
r1.YData=15:35;

subplot(3,3,7)
plot(flip(age),full_series,'Color','black');
hold on
plot(flip(age(1:(n_months*215))),flip(bio_series3(1:(n_months*215),1)),'Color','0.8 0.8 0.8');
hold on
x2 = [min(min(age_res(:,201:215))) max(max(age_res(:,201:215))) max(max(age_res(:,201:215))) min(min(age_res(:,201:215)))];
y2 = [min(min(full_series3_res(:,201:210)))    min(min(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))];
patch(x2,y2,'red','EdgeColor','red')
hold on
x2 = [min(min(age_res(:,201:210))) max(max(age_res(:,201:210))) max(max(age_res(:,201:210))) min(min(age_res(:,201:210)))];
y2 = [min(min(full_series3_res(:,201:210)))    min(min(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))  max(max(full_series3_res(:,201:210)))];
patch(x2,y2,[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250],'FaceColor','none')

xlabel('Months');
ylabel('Temperature (°C)')
ylim([15 35])

subplot(3,3,8)

for i=1:1:length(v3a) % after bio 
    plot(v3a(i,:),p3a(i,:),'Color','[1, 0, 0, 0.03]')
    hold on
end

hold on

for i=1:1:length(v3b) % before bio
    plot(v3b(i,:),p3b(i,:),'Color','[0.9290 0.6940 0.1250 0.03]')
    hold on
end

hold on
%plot(mean(v3b),mean(p3b),'Color','[0.9290 0.6940 0.1250]','LineWidth',2) % mean kernels before bio
hold on
%plot(mean(v3a),mean(p3a),'Color','red','LineWidth',2) % mean kernels after bio

xlabel('Temperature (ºC)');
ylabel('Prob.Density');
xlim([15 35])
ylim([0 0.7])

subplot(3,3,9)
errorbar(median(qx3,2),median(qy3,2),serrY3(1,:),serrY3(2,:),'linestyle','none','color',[0.9290 0.6940 0.1250],'capsize',0,'linewidth',1.5);
hold on
errorbar(median(qx3,2),median(qy3,2),serrX3(1,:),serrX3(2,:),'horiz','color',[0.9290 0.6940 0.1250],'linestyle','none','capsize',0,'linewidth',1.5);
hold on
plot(median(qx3,2),median(qy3,2),'kd-','MarkerSize',6,'MarkerFaceColor','k','LineStyle','none');
xlabel('Non-Bioturbated (ºC)','Color','[0.9290 0.6940 0.1250]');
ylabel('After Bio (ºC)','Color','r');
rl = refline(1,0); 
rl.Color = 'k';
rl.LineWidth = 1.5;
set(gca, 'YAxisLocation', 'right')
xlim([21.5 29.5])
ylim([21.5 29.5])
r1.XData=21:31;
r1.YData=21:31;

annotation('textbox', [0.08,0.77,0.21,0.21], 'String', 'a)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.77,0.21,0.21], 'String', 'b)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.77,0.21,0.21], 'String', 'c)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.47,0.21,0.21], 'String', 'd)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.47,0.21,0.21], 'String', 'e)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.47,0.21,0.21], 'String', 'f)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.08,0.19,0.21,0.21], 'String', 'g)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.36,0.19,0.21,0.21], 'String', 'h)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('textbox', [0.64,0.19,0.21,0.21], 'String', 'i)', 'FontSize', 9, 'EdgeColor', 'none');
annotation('arrow',[0.32 0.32],[0.72 0.76],'Color','blue','LineWidth',1)
annotation('arrow',[0.32 0.32],[0.93 0.89],'Color','red','LineWidth',1)
annotation('textarrow',[0.32 0.28],[0.69 0.69],'String','Bioturbation')
annotation('arrow',[0.255 0.255],[0.42 0.46],'Color','[0.4660 0.6740 0.1880]','LineWidth',1)
annotation('arrow',[0.255 0.255],[0.65 0.61],'Color','red','LineWidth',1)
annotation('arrow',[0.19 0.19],[0.13 0.17],'Color','[0.9290 0.6940 0.1250]','LineWidth',1)
annotation('arrow',[0.19 0.19],[0.34 0.30],'Color','red','LineWidth',1)

clear;
