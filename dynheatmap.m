function results=dynheatmap(seas,IFAtop_mage,IFApast_mage,run,ageS,m_zS,BDin,abuIFA,abu,ageAbu,le,clim_startm,scale)
%                          
% Function that bioturbates a longer time series in order to calculate the
% probability of pick individuals older than IFA time slice in dynamical 
% scenarios (i.e., dynamical SAR, SMLD and or abundance).
%
% --> The mixing model builds on the SEAMUS software bioturbation 
% 
% Source SEAMUS:
% https://github.com/bryanlougheed/seamus
% Lougheed, B. C. (2020). Seamus (v1. 20): a δ14 c-enabled, single-specimen sediment
% accumulation simulator. Geoscientific Model Development, 13 (1), 155–168
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 20-Dec-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

le=le*1000; % in years
Trace_full=le*12 -1 + clim_startm:-1:clim_startm; % last le ka monthly time series, relative to 1950CE

if ageS(end)<(length(Trace_full)- clim_startm/12) % add cte abu ageS(end) to lenght(Trace_full)
        m_zSe=(m_zS(length(m_zS))-m_zS(length(m_zS)-1))*((length(Trace_full)/12)-ageS(length(ageS)))...
        /(ageS(length(ageS))-ageS(length(ageS)-1));
        m_zS=[m_zS,m_zS(length(m_zS)) + m_zSe];
        ageS=[ageS,(length(Trace_full)- clim_startm)/12];
        abuIFA=[abuIFA,abuIFA(length(abuIFA))];
        BDin=[BDin,BDin(length(BDin))]; 
        if abu==2
            ageAbu=[ageAbu,(length(Trace_full)- clim_startm)/12];
        end
end

  
% get depths idx per month
cum_month=(clim_startm:1:length(Trace_full)-1 + clim_startm)/12; 
depth_ind=(interp1(ageS,m_zS, cum_month,'linear'));
depth_ind(depth_ind <0)=depth_ind(depth_ind <0) -1; % correct -0cm vs +0cm assign diff depths
depths_ind=flip(depth_ind);
pos_depth=fix(depths_ind);  % assign same depth integerer to each discrete depth
z=max(pos_depth):-1:min(pos_depth);  % final core depths 
BD=interp1(m_zS,BDin,z,'linear','extrap'); % dowcore BD (cm)
age=(length(Trace_full)-1 + clim_startm:-1:clim_startm);

if range(abuIFA) ~= 0  % dynamic abundance
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
        seas_indx=1:1:length(Trace_full);
        ll_res=reshape(seas_indx,12,[]);
        ll_seas=reshape(ll_res(seas,:),[],1); % index seas months
        rep_ind=rep_ind(ll_seas);
        
        % take seas data and repeat it according to abu
        depths_ind=repelem(depths_ind(ll_seas),rep_ind);
        age=repelem(age(ll_seas),rep_ind);
        pos_depth=repelem(pos_depth(ll_seas),rep_ind);
            
    else % yearly specimens
         depths_ind=repelem(depths_ind,rep_ind);
         age=repelem(length(Trace_full):-1:1,rep_ind); 
         pos_depth=repelem(pos_depth,rep_ind);
    end

else  % no dynamical abundance
    
    if length(seas)<12  % if seasonal preference  
         % Find indx for seas months (before bio). Note that simulation ends at
         % december (Present) 
         seas_indx=1:1:length(Trace_full);
         ll_res=reshape(seas_indx,12,[]);
         ll_seas=reshape(ll_res(seas,:),[],1); %index seas months
         depths_ind=depths_ind(ll_seas);
         age=age(ll_seas);
         pos_depth=pos_depth(ll_seas);
            
    end
end

% mixing model
age_bio=repmat(age',1,run);
z=flip(unique(pos_depth));

if any(z<0)
    z=[z,z(end)-1];
    BD=[BD, BD(end)];
end

for r= 1:run
    for i = 1:1:length(z)
        ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD(i)); % find indx of depths values inside
        rndind = ind(randperm(length(ind))); % random displacement
        age_bio(ind,r)=age_bio(rndind,r); % bioturbated age time series
    end
end

% Get IFA time slices data
     
% order from younger to older;
age_bio=flip(age_bio);
pos_depth=flip(pos_depth);
z=flip(z);

if any(z<0)
    z=z(2:end);
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

% Bioturbated CoreTop data
age_bio_top=cell(length(IFAtop_mage),run);

for ii=1:length(IFAtop_mage)
    [~,idx]=min(abs(med_bio-IFAtop_mage(ii))); 
    for i=1:run
        age_bio_top{ii,i} = age_bio(s_idx(idx(i)):e_idx(idx(i)),i); 
    end
end

% Bioturbated downcore data
age_bio_Past=cell(length(IFApast_mage),run); % bioturbated ages downcore IFA data depths

for ii=1:length(IFApast_mage)
    [~,idx]=min(abs(med_bio-IFApast_mage(ii))); 
    for i=1:run
        age_bio_Past{ii,i} = age_bio(s_idx(idx(i)):e_idx(idx(i)),i); 
    end
end

results.age_bioheat_Past=age_bio_Past; % bioturbated ages IFA downcore 
results.age_bioheat_top=age_bio_top;   % bioturbated ages IFA coretop
end                                    
