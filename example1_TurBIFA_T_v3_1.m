%% Example 1. TurBIFA model

%% About
%
% This example uses Mg/Ca-T IFA data (in ºC) and constant SAR,BD and abundance
% 
% Note: just run the script to explore TurBIFA_T outputs
%
% --> TurbIFA_T mixing model builds on the SEAMUS software bioturbation 
% --> TurbIFA_T IFA builds on the INFAUNAL and QUANTIFA algorithms
%     to estimate the noise in the std and shape (q-q space) 
%     of the IFA distribution, respectively.
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
% Source INFAUNAL
% https://www.mathworks.com/matlabcentral/fileexchange/47695-infaunal
% https://doi.org/10.1002/palo.20037
%
% Grafical Outputs:
% *************************************************************************
% 1. A heatmap chart, which displays the probability of picking in 
% a given discrete depth sample 1 to 20 individuals that precipitated 
% their shells during time slices 1-10 kyr older than the 
% IFA sample time slice.
%
% 2. Comparisons of core-top versus downcore IFA distributions in the 
% quantile-quantile space with 95% confidence intervals (based on QUANTIFA
% software matlab code for q-q graphical output)
%.
% 3. Summary table of the uncertainties of the IFA Mg/Ca distribution, 
% and of the probability that the change in the standard deviation (IFA 1σ) of 
% IFA distribution is significant
%
% 4. Probability significant increase / decrease variability for
% pseudo IFA datasets resampled in a non-bioturbated vs bioturbated sediment 
% archive (qq space)
%
% Other m-files required: quantifaerr,analiBH
% 
% MAT-files required: - data_examplesTurBIFA.mat
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 27-April-2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Input data (run as it is)

Trace = 1;              % bioturbate Trace data (SAT)
run = 1;                % number of bioturbated synthetic sediment cores (default is 50, now 1 just to run faster)
abu = 0;                % constant abundance
anerr = 0.5;            % analytical uncertainty (ºC)
exerr = 0;              % extra gaussian noise from e.g., calibration (1σ,ºC) 
IFAtop_mage = 1750;     % mean age Coretop IFA (years, BP)
ystart_top = 2000;      % start year coretop IFA data timeslice (years, BP)
yend_top = 1500;        % end year  coretop IFA data timeslice (years, BP)
IFApast_mage = 5750;    % mean age downcore IFA (years, BP)
ystart_past = 6000;     % start year past IFA data time-slice (years, BP)
yend_past = 5500;       % end year past IFA data time-slice (years, BP)
dynamic = 0;            % constant SAR,BD and abundance
seas = 1:12;            % yearly specimens
sed_r = 2;              % sedimentation rate (cm/kyr)
BD = 5;                 % bioturbation depth (cm)
num_q = 50;             % numer of quantiles;
mc = 1000;              % repeat "picking" mc times
clim_start = -40;       % relative to 1950 CE
months_kyrs = 12*1000;

%% 1.1 Load SAT TraCE data  (lat= -0.13ºN, lon= 204E)

Trace_full=load('data_examplesTurBIFA.mat','TraceSAT_ex1');
Trace_full=Trace_full.TraceSAT_ex1;

%% 1.2 Load IFA data 

X=load('data_examplesTurBIFA.mat','LH_data1');
X=X.LH_data1;
Y=load('data_examplesTurBIFA.mat','Past_data1');
Y=Y.Past_data1;

% IFA time slices in months
clim_startm= (clim_start*12);  % start month climate time series, relative to 1950CE

%% 2 Prepare TºC time series before bioturbation simulation

if dynamic==0

    %% 2.1 Constant SAR,SMLD & individuals abundance

    z=flip(1:1:(length(Trace_full)/months_kyrs)*sed_r); % depths                              
    n_months=round(months_kyrs/sed_r);       % nº months (have to be round) at each discreate depth
    pos_depth=repelem(z,n_months); 
    depths_ind=(pos_depth(1)+1)-1/n_months:-1/n_months:(z(end));
    
    if length(pos_depth)> length(Trace_full)   % to avoid error due to round z & n_months
          ln=length(pos_depth);
          lt=length(Trace_full);
          depths_ind=depths_ind(((ln-lt)+1):end);
          pos_depth=pos_depth(((ln-lt)+1):end);
    elseif length(pos_depth)< length(Trace_full) % to ensure mixing most recent months
          ln=length(pos_depth);
          lt=length(Trace_full);
          depths_ind=[repelem(depths_ind(1),(lt-ln)),depths_ind];
          pos_depth=[repelem(z(1),(lt-ln)),pos_depth]; 
    end
    
    BD=repelem(BD,length(z),1);
    age=(length(Trace_full)-1 + clim_startm:-1:clim_startm);

    if length(seas)<12   % if seasonal preference
         % Find indx for seas months (before bio). Note that simulation ends at
         %  december (Present)
         seas_indx=1:length(Trace_full);
         ll_res=reshape(seas_indx,12,[]);
         ll_seas=reshape(ll_res(seas,:),[],1); % index seas months
         Trace_full=Trace_full(ll_seas);
         depths_ind=depths_ind(ll_seas);
         pos_depth=pos_depth(ll_seas);
         age=age(ll_seas);
    end

else  
     %% 2.2 Dynamical SAR and/or SMLD and/or IFA abundance

     if ageS(end)<(length(Trace_full)/12) - clim_startm % add cte abu ageS(end) to lenght(Trace_full)
            m_zSe=(m_zS(length(m_zS))-m_zS(length(m_zS)-1))*((length(Trace_full)/12)-ageS(length(ageS)))...
            /(ageS(length(ageS))-ageS(length(ageS)-1));
            m_zS=[m_zS,m_zS(length(m_zS)) + m_zSe];
            ageS=[ageS,(length(Trace_full)- clim_startm)/12];
            BDin=[BDin,BDin(length(BDin))];   
            if abu==1 | abu==0
                abuIFA=[abuIFA,abuIFA(length(abuIFA))];
            end
     end

     if ageS(1)>clim_start  % add cte abu from ageS(1) until core top
             m_zSe=-1*(m_zS(2)-m_zS(1))*(ageS(1)-(clim_start/1000))...
            /(ageS(2)-ageS(1));
             m_zS=[m_zSe,m_zS];
             ageS=[clim_start,ageS];
             BDin=[BDin(1),BDin];
             if abu==1 | abu==0
                 abuIFA=[abuIFA(1),abuIFA];
             end
     end

     if abu==2
         if ageAbu(end)<(length(Trace_full)/12) % add cte abu ageAbu(end) to lenght(Trace_full)
             ageAbu=[ageAbu,(length(Trace_full)- clim_startm)/12];
             abuIFA=[abuIFA,abuIFA(length(abuIFA))];  
           
         end
         if ageAbu(1)>clim_start  % add cte abu from ageAbu(1) until core top
             ageAbu=[clim_start,ageAbu];
             abuIFA=[abuIFA(1),abuIFA];
             
         end
     end
  
     % get depths idx per month
    
     cum_month= (clim_startm:1:length(Trace_full)-1 + clim_startm)/12;
     depth_ind=(interp1(ageS,m_zS, cum_month,'linear','extrap'));
     depth_ind(depth_ind <0)=depth_ind(depth_ind <0) -1; % correct -0cm vs +0cm assign diff depths
     depths_ind=flip(depth_ind);
     pos_depth=fix(depths_ind);  % assign same depth integerer to each discrete depth
     z=max(pos_depth):-1:min(pos_depth);  % final core depths 
          
     BD=interp1(m_zS,BDin,z,'linear','extrap'); % downcore BD (cm)
     age=(length(Trace_full)-1 + clim_startm:-1:clim_startm);  

     if range(abuIFA) ~= 0  % if abuIFA is not constant
         % scale abundance 
         max_abu=max(abuIFA);
         sabuIFA=round((abuIFA*scale)/max_abu); % scale to repeat monthly values max scale times
         if abu==1
             rep_ind=flip(round(interp1(ageS,sabuIFA,cum_month,'linear','extrap')));
         else % abu==2
             rep_ind=flip(round(interp1(ageAbu,sabuIFA,cum_month,'linear','extrap')));
         end
         rep_ind(rep_ind<0)=0;

         if length(seas)<12  % if seasonal preference
             % Find indx for seas months (before bio). Note that simulation ends at
             % december (Present)
             seas_indx=1:1:length(age);
             ll_res=reshape(seas_indx,12,[]);
             ll_seas=reshape(ll_res(seas,:),[],1); % index seas months
             rep_ind=rep_ind(ll_seas);
            
             % take seas data and repeat it according to abu
             Trace_full=repelem(Trace_full(ll_seas),rep_ind);
             depths_ind=repelem(depths_ind(ll_seas),rep_ind);
             age=repelem(age(ll_seas),rep_ind);
             pos_depth=repelem(pos_depth(ll_seas),rep_ind);
            
         else  % yearly specimens
             depths_ind=repelem(depths_ind,rep_ind);
             age=repelem(age,rep_ind);
             Trace_full=repelem(Trace_full,rep_ind);
             pos_depth=repelem(pos_depth,rep_ind);
         end

     else  % no dynamical abundance

         if length(seas)<12     % if seasonal preference  
             % Find indx for seas months (before bio). Note that simulation ends at
             % december (Present) 
             seas_indx=1:1:length(Trace_full);
             ll_res=reshape(seas_indx,12,[]);
             ll_seas=reshape(ll_res(seas,:),[],1); %index seas months
             Trace_full=Trace_full(ll_seas);
             depths_ind=depths_ind(ll_seas);
             age=age(ll_seas);
             pos_depth=pos_depth(ll_seas);
            
         end
     end
end

%% 3. Sediment Mixing Model

z=flip(unique(pos_depth));  % final core depths

if any(z<0) % ensure mixing coretop
    z=[z,z(end)-1];
    BD=[BD, BD(end)];
end

age_bio=repmat(age',1,run); 
Trace_bio=repmat(Trace_full,1,run);

for r= 1:run
    for i = 1:length(z)
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD(i)); % find indx of depths values inside
        rndind = ind(randperm(length(ind))); % random displacement
        Trace_bio(ind,r)=Trace_bio(rndind,r); % bioturbated TºC time series
        age_bio(ind,r)=age_bio(rndind,r); % bioturbated age time series
    end
end

%% 3.1 Get IFA time slices data
     
% order from younger to older
Trace_bio=flip(Trace_bio); 
Trace_full=flip(Trace_full);
age=flip(age);
age_bio=flip(age_bio);
z=flip(z);
pos_depth=flip(pos_depth);

% IFA time slices in months
smonth_top=ystart_top*12 + clim_start*12;       % start month of LH (or coretop) IFA data
emonth_top= (yend_top*12)-11 + clim_start*12;   % end month of LH (or coretop) IFA  data

smonth_past=ystart_past*12 + clim_start*12;     % start month of downcore IFA data
emonth_past= (yend_past*12)-11 + clim_start*12; % end month of downcore IFA data

if any(z<0)
    z=z(2:end);  % not needed, just for bioturbation simulation
end

% Calculate mean age for each discrete depth inside the bioturbated synthetic
% core sediment archive.
s_idx=zeros(length(z),1); % first indx for each discrete depth
e_idx=zeros(length(z),1); % last indx for each discrete depth
med_bio=zeros(length(z),run); % mean age per cm after bio

for i = 1:length(z)
    pos = find(pos_depth == z(i));
    s_idx(i)=min(pos); 
    e_idx(i)=max(pos); 
    med_bio(i,:)=mean(age_bio(s_idx(i):e_idx(i),:))/12; 
end

% Non-bioturbated CoreTop data
[~,idx_emonthtop]=min(abs(age-emonth_top)); % index (position) emonth_top 
[~,idx_smonthtop]=min(abs(age-smonth_top)); % index smonth_top 
top_nonbio=Trace_full(idx_emonthtop:idx_smonthtop); % CoreTop IFA time slice non bioturbated ("True variability") 
  
% Non-bioturbated downcore data 
[~,idx_emonthpast]=min(abs(age-emonth_past)); % index emonth_past      
[~,idx_smonthpast]=min(abs(age-smonth_past)); % index smonth_past
Past_nonbio=Trace_full(idx_emonthpast:idx_smonthpast); % Downcore IFA time slice non bioturbated

% Bioturbated CoreTop data
top_bio=cell(length(IFAtop_mage),run); % bioturbated coretop IFA data

for ii=1:length(IFAtop_mage)
    [~,idx]=min(abs(med_bio-IFAtop_mage(ii))); 
    for i=1:run
        top_bio{ii,i} = Trace_bio(s_idx(idx(i)):e_idx(idx(i)),i); 
    end
end

% Bioturbated downcore data
Past_bio=cell(length(IFApast_mage),run); % bioturbated downcore IFA data 

for ii=1:length(IFApast_mage)
    [~,idx]=min(abs(med_bio-IFApast_mage(ii))); 
    for i=1:run
        Past_bio{ii,i} = Trace_bio(s_idx(idx(i)):e_idx(idx(i)),i); 
    end
end
 

%% 4 "Pick" data and calculate noise 

% 1. Non-bioturbated time slice
pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,top_nonbio,Past_nonbio,run);

errTnbio = pick_output.Xqq_err;     % uncertainty non-bio Core top population (qq space)
errPnbio= pick_output.Yqq_err;      % uncertainty non-bio downcore population (qq space)
qTnbio=pick_output.qqx;             % core top non-bio resampled dataset (qq-space, mean removed) 
qPnbio=pick_output.qqy;             % downcore non-bio resampled dataset (qq-space, mean removed) 

% 2. Bioturbated time slice
pick_output=quantifaerrv3(X,Y,anerr,exerr,num_q,top_bio,Past_bio,run,IFAtop_mage,IFApast_mage);

errTbio = pick_output.Xqq_err;   % uncertainty After bio Core top population (qq space)
errPbio= pick_output.Yqq_err;    % uncertainty After bio Downcore population (qq space)
sigmaX_per =pick_output.Xstd_per;% uncertainty Core top IFA 1-sigma distribution value (95% CI)
sigmaY_per =pick_output.Ystd_per;% uncertainty downcore IFA 1-sigma distribution value (95% CI)
std_i=pick_output.std_i;         % probability significant increase Y population 1-sigma value
std_d=pick_output.std_d;         % probability significant decrease Y population 1-sigma value
qTbio=pick_output.qqx;           % core top bio resampled dataset (qq-space, mean removed) 
qPbio=pick_output.qqy;           % downcore bio resampled dataset (qq-space, mean removed) 


%% 4.1 Probability of pick n individuales older than x

if dynamic==0 && length(X) == length(Y) % equal probability along core if SAR, SMLD and nº IFA forams cte
   
    % calculate heatmap with synthetic age distribution (analytical B & H solution)
    data=analiBH(sed_r,BD,X);
    
    % probabality of pick 1-20 indiv of 1-10 kyrs older than (mean) age 1cm
    pbt=data.pp;
       
else
    %  to avoid biases due to high BD/SAR
       
    if dynamic==0  % length(X) ~= length(Y)
        
         % X population
         data=analiBH(sed_r,BD,X);
         tpbt=data.pp;

         % Y population
         data=analiBH(sed_r,BD,Y);
         ppbt=data.pp;
        
    else
        results=dynheatmap(seas,IFAtop_mage,IFApast_mage,run,ageS,m_zS,BDin,abuIFA,abu,ageAbu,le,clim_startm,scale);
        
        age_bio_Pastolka= results.age_bioheat_Past;
        age_bio_topolka = results.age_bioheat_top;
        num_bio=1:20;     % number of bioturbated individuals 
        age_offset=10;    % age offset (kyr)
        top=zeros(mc,age_offset,length(num_bio),run);
        past=top;
        num_bio_ex=reshape(repelem(num_bio,mc,age_offset),mc,age_offset,[]);
        ageoffset_eX=reshape(repelem(1:age_offset,length(X),mc),length(X),mc,[]);
        ageoffset_eY=reshape(repelem(1:age_offset,length(Y),mc),length(Y),mc,[]);
       
       for r=1:run
           %core top
           if length(IFAtop_mage) ==1 % IFA data picked 1 depth
               mage_coretop = mean(cell2mat(age_bio_topolka(r)));
               cage_bio_top = cell2mat(age_bio_topolka(r));
           else % more than 1 depth
               mage_coretop=mean(reshape(vertcat(age_bio_topolka{1:length(IFAtop_mage),r}),[],1));
               cage_bio_top = reshape(vertcat(age_bio_topolka{1:length(IFAtop_mage),r}),[],1);
           end
           % picking 
           idx = randi(length(cage_bio_top),length(X),mc);
           older_ind= cage_bio_top(idx)> (ageoffset_eX*months_kyrs) + mage_coretop;
           n_older_ind= squeeze(sum(older_ind)) > num_bio_ex;
           top(:,:,:,r)= n_older_ind;

           % downcore slice
           if length(IFApast_mage) ==1 
               mage_past = mean(cell2mat(age_bio_Pastolka(r)));
               cage_bio_Past = cell2mat(age_bio_Pastolka(r));
           else
               mage_past=mean(reshape(vertcat(age_bio_Pastolka{1:length(IFApast_mage),r}),[],1));
               cage_bio_Past = reshape(vertcat(age_bio_Pastolka{1:length(IFApast_mage),r}),[],1);
           end
           % picking 
           idx = randi(length(cage_bio_Past),length(Y),mc);
           older_ind= cage_bio_Past(idx)> (ageoffset_eY*months_kyrs) + mage_past;
           n_older_ind= squeeze(sum(older_ind)) > num_bio_ex;
           past(:,:,:,r)=n_older_ind;
       end
    
       % probabality of pick 1-20 indiv of 1-10 kyrs older than (mean) Top IFA
       % time slice
       tpbt=100*squeeze(sum(sum(top,1),4)/(mc*run));  % coretop
       ppbt=100*squeeze(sum(sum(past,1),4)/(mc*run)); % downcore 
   
    end
end

%% 4.2 Probability significant increase / decrease variability for a dataset resampled in a
% non-bioturbated vs bioturbated sediment archive (qq space)

% non bioturbated
diff_nbio=abs(qPnbio)-abs(permute(qTnbio, [1, 3, 2]));
% significant quantiles above 1:1 line
qq_ab=qPnbio>(permute(qTnbio, [1, 3, 2])) & abs(diff_nbio) > abs(errPnbio(1,:)') & abs(diff_nbio) > abs(errTnbio(2,:)');
% significant quantiles behind 1:1 line
qq_be=qPnbio<(permute(qTnbio, [1, 3, 2])) & abs(diff_nbio) > abs(errPnbio(2,:)') & abs(diff_nbio) > abs(errTnbio(1,:)');
%increase
incres = qq_ab & diff_nbio > 0 | qq_be & diff_nbio >0;
%decrease
decres = qq_be >0 & diff_nbio < 0 | qq_ab >0 & diff_nbio <0;
% prob significant increase & decrease (non-bioturbated)
signbio_in=100*(sum(sum(incres,3),2))/(mc*mc);  % prob increase nbio
signbio_de=100*(sum(sum(decres,3),2))/(mc*mc);  % prob decrease nbio

% Declare data (Bioturbated)
sigbio_in=zeros(num_q,run);   % prob increase bio
sigbio_de=zeros(num_q,run);   % prob decrease bio
aberrPbio=abs(errPbio);
aberrTbio=abs(errTbio);

for r=1:run
    qTbior=qTbio(:,:,r);
    qPbior=qPbio(:,:,r);
    diff_bio=abs(qPbior)-abs(permute(qTbior, [1, 3, 2])); %difference
    abdiff_bio=abs(diff_bio);
    % significant quantiles above 1:1 line
    qq_ab=qPbior>(permute(qTbior, [1, 3, 2])) & (abdiff_bio) > (aberrPbio(1,:,r)') & (abdiff_bio) > (aberrTbio(2,:,r)'); 
    % significant quantiles behind 1:1 line
    qq_be=qPbior<(permute(qTbior, [1, 3, 2])) & (abdiff_bio) > (aberrPbio(2,:,r)') & (abdiff_bio) > (aberrTbio(1,:,r)'); 
    %increase
    incres = qq_ab & diff_bio > 0 | qq_be >0 & diff_bio >0 ;
    %decrease
    decres = qq_be >0 & diff_bio < 0 | qq_ab >0 & diff_bio <0;
    % prob significant increase & decrease (bioturbated)
    sigbio_in(:,r)=100*(sum(sum(incres,3),2))/(mc*mc);  % prob increase bio
    sigbio_de(:,r)=100*(sum(sum(decres,3),2))/(mc*mc);  % prob decrease bio 
end
        
%% 5. TurBIFA Outputs

% Output 1. Heatmap chart, which displays the probability of picking in a given 
% depth sample 1 to 20 individuals that precipitated their shells during 
% time slices 1-to 10 kyr older than the IFA sample time slice

figure

if dynamic==0 && length(X) == length(Y)
    h1=heatmap(1:1:20,flip(1:10),flip(pbt));
    h1.Colormap=hot;
    xlabel('nº individuals')
    ylabel('Age difference (kyr)')

else
    subplot(1,2,1)
    h1=heatmap(1:1:20,flip(1:10),flip(tpbt));
    h1.Colormap=hot;
    xlabel('nº individuals')
    ylabel('Age difference (kyr)')
    title('Coretop IFA')

    subplot(1,2,2)
    h2=heatmap(1:1:20,flip(1:10),flip(ppbt));
    h2.Colormap=hot;
    xlabel('nº individuals')
    ylabel('Age difference (kyr)')
    title('Downcore IFA')
end

% Output 2. Core-top versus Downcore IFA distributions in the quantile-quantile
% space with 95% confidence intervals 

qx=quantile(X-mean(X),num_q); % Coretop mean removed qq
qy=quantile(Y-mean(Y),num_q); % Downcore mean removed qq

figure
errorbar(qx,qy,median(errPbio(1,:,:),3),median(errPbio(2,:,:),3),...
    median(errTbio(1,:,:),3),median(errTbio(2,:,:),3),...
    'o','Color','r','MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',0);     
rl = refline(1,0); 
rl.Color = 'k';
xlabel('CoreTop IFA  Quantiles  (ºC, mean removed)');
ylabel('Downcore IFA Quantiles (ºC, mean removed)');

                      
% Output 3. Summary table of uncertainties (1σ-IFA) and of the probability that changes
% in the standard deviation (1σ) of IFA distribution are significant
                   
figure
LastName = {'IFAσ';'Error (P2.5)';'Error (P97.5)'};
Coretop = [std(X);median(sigmaX_per(1,:));median(sigmaX_per(2,:))];
Past = [std(Y);median(sigmaY_per(1,:));median(sigmaY_per(2,:))];
T = table(Coretop,Past,'RowNames',LastName);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
hold on
LastName = {'pbb sig increase (%)'; 'pbb sig decrease (%)'};
prob_sig=[median(std_i); median(std_d)];
TT=table(prob_sig,'RowNames',LastName);
uitable('Data',TT{:,:},'ColumnName',TT.Properties.VariableNames,...
    'RowName',TT.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 0.7]);

% Output 4. Probability significant increase / decrease variability for
% pseudo IFA datasets resampled in a non-bioturbated vs bioturbated sediment 
% archive (qq space)

figure 
subplot(1,3,1)
bio=scatter(qTbio(:,:,:,1),qPbio(:,:,:,1),'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.4);
hold on
nobio=scatter(qTnbio,qPnbio,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.4);
rl = refline(1,0); 
rl.Color = 'k';
xlabel('Coretop (ºC, mean removed)')
ylabel('Downcore (ºC, mean removed)')
legend([bio(1),nobio(2)],'resampled After Bio','resampled Before Bio')

% increase
subplot (1,3,2)
in=[median(signbio_in,2)'; median(sigbio_in,2)'];
bar(in','EdgeColor','none')
xlabel('quantile')
ylabel('probability (%)')
ylim([0 100])
legend('Before Bio','After Bio')
title('Downcore increase variability')

subplot(1,3,3)
de=[median(signbio_de,2)'; median(sigbio_de,2)'];
bar(de','EdgeColor','none')
xlabel('quantile')
ylabel('probability (%)')
ylim([0 100])
legend('Before Bio','After Bio')
title('Downcore decrease variability')

%% 6. Save data outputs

outputs.qx=qx;                     
outputs.qy=qy;                     
outputs.errTnbio=errTnbio;          
outputs.errPnbio=errPnbio;         
outputs.errTbio=errTbio;            
outputs.errPbio=errPbio;            
outputs.errsigmaX=sigmaX_per;       
outputs.errsigmaY=sigmaY_per;       
outputs.std_i=std_i;                
outputs.std_d=std_d;                
outputs.qTnbio=qTnbio;              
outputs.qPnbio=qPnbio;              
outputs.qTbio=squeeze(qTbio);     
outputs.qPbio=squeeze(qPbio);       
outputs.qqsignbio_in=signbio_in;    
outputs.qqsignbio_de=signbio_de;    
outputs.qqsigbio_in=sigbio_in;      
outputs.qqsigbio_de=sigbio_de;     
if dynamic==0 && length(X) == length(Y)
    outputs.pbt = pbt;
else
    outputs.tpbt=tpbt;              
    outputs.ppbt=ppbt;              
end

clearvars -except outputs % clear all vars except outputs


 