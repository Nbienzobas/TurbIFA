% 
%% About
%
% bioprob is a matlab script that calculate the probability of pick individuales that grew
% during time slices older than the mean age of the IFA data depth slice
%
% ---> if the sed.rate, bioturbation depth and abundance are constant,
%      bioprob calculates the exponential age distribution using an scale
%      parameter equal to BD/SAR
%
%      Inputs
%      *******************************************************************
%      - size of the IFA dataset, that is, number of IFA individuals
%       meausured
%      - mean sedimentation rate (cm/kyr) 
%      - mean bioturbation depth (cm)
%      
% ---> if one or more parameters varies over time, bioprob perform a random
%      step by step bioturbation simulation based on SEAMUS software.
%
%      Source SEAMUS:
%      https://github.com/bryanlougheed/seamus
%      Lougheed, B. C. (2020). Seamus (v1. 20): a δ14 c-enabled, single-specimen sediment
%      accumulation simulator. Geoscientific Model Development, 13 (1), 155–168
% 
%      Inputs for dynamical scenarios
%      *******************************************************************
%     - size of the IFA dataset
%     - constant sedimentation rate (cm/kyr) or dynamical sedimentation rate
%      (insert data as age-discrete rounded depth)
%     - constant bioturbation depth (cm) or dynamical bioturbation depth (cm)
%     - abu= 0  Constant foraminiferal abundance through time (default)
%     - abu = 1 downcore changes in foraminiferal abundance (per depth). Abundance data derived 
%                from individual census counts. Data will be scaled
%                to an user defined scale parameter e.g., 10. As the scale paramter increases, smaller changes
%                in the IFA abundance can be simulated) 
%     - abu = 2 Temporal changes foraminiferal abundance (e.g., relative flux. Data will be scaled
%               to an user defined scale parameter e.g., 10. As the scale paramter increases, smaller changes
%               in the IFA abundance can be simulated)  
%     - run = number of bioturbated synthehic cores (50 by default)**
%     - mean age Past IFA depth = e.g. 311 ka. NOTE: If more than 1 discrete depth
%                                 picked for IFA, insert the mean age of each one (e.g, 311 311.25 311.5 ...)
%     - Simulation start year = e.g. 400 ka
%     - Simulation end year   = e.g. 300 ka
%     - individuals seasonality = e.g., yearly ('Y'); summer season (6 7 8)...
%
%     * Note that the IFA time slice should be closer to the end of the simulation 
%      (younger age) than to the beginning of the simulation (oldest age).
%     ** Note: For a scale parameter value >10, I recommend to decrease the number of bioturbated
%      synthehic sediment core archived ("run" variable) to avoid memory issues.
% 
% Output
% *************************************************************************
% pp= probability matrix size 10x20 where columns = age_offset and 
% rows = number of individuals within age_offset range
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 27-April-2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. Constant SAR, BD and/or abundance?

dynamic = questdlg('Constant sedimentation rate, Bioturbation depth and IFA specie abundance?', ...
	'Dynamical data', ...
	'Yes','No','Yes');

if ismember('Yes',dynamic)==1
    dynamic=0;  % constant SAR, BD and abu
else
    dynamic=1;  % dynamic SAR,BD and/or abu
end

if dynamic==0
    
    %% 1.1a. Input constant parameters
    x = inputdlg({'mean sedrate cm/kyrs:','Bioturbation depth (cm):','nº foraminfera:'});
    sed_r=str2double(x{1,1});          % mean sed.rate at core location
    BD=str2double(x{2,1});             % Bioturbation depth (cm)
    num_p=str2double(x{3,1});          % number of foraminifera picked for IFA


    %% 1.2a. Calculate probability
    % calculate heatmap with synthehic age distribution (analytical B & H solution)
    indivs=repelem(1,num_p,1);
    data=analiBH(sed_r,BD,indivs);
    
    % probabality of pick 1-20 indiv of 1-10 kyrs older than (mean) age 1cm
    pp=data.pp;
           

else  % dynamical SAR , BD and/or abundance

    %% 1.1b. Input dynamical parameter(s)
    x = inputdlg({'nº foraminfera:','Simulation start year (e.g., 400 ka):', ...
    'Simulation end year (e.g., 300 ka):','mean IFA age (e.g.,311 ka):'...
    ,'abu (e.g., 0)','run'},...
              'Model Parameters');

    num_p =str2double(x{1,1});                  % size IFA dataset
    sstart_y=str2double(x{2,1})*1000;           % simulation start year e.g., 400 ka
    send_y=str2double(x{3,1})*1000;             % simulation end year e.g., 300ka
    IFApast_mage=sscanf(x{4,1},'%f')*1000;      % mean age IFA discrete depth
    abu=str2double(x{5,1});                     % abundance simulation      
    run=str2double(x{6,1});                     % number of bioturbated synthehic cores

    age=((sstart_y-1)*12:-1:(send_y*12)-11);    % monthly time series
     
    % Seasonal Preference?
    prompt = {'Year= Y;  DJF= 1 2 12; JJA= 6 7 8;'};
    dlgtitle = 'Seasonal Preference?';
    answer = inputdlg(prompt,dlgtitle,[1, 50]);
    
    if ismember('Y',answer)==1 % to avoid the use of stru2num
        seas=1:12;
    else
    seas=sscanf(answer{1, 1},'%f')'; 
    end

    if isempty(seas)
        error('The variable "seas" is empty. This could be due to the assigment of a non numerical value')
    end

    if abu == 1 || abu==0 % downcore foraminiferal individual census counts or cte
        
        prompt = {'Enter core depths (in cm !!):','Enter Ages (kaBP):','Enter Bioturbation depth (cm):', 'Enter IFA specie abudance'};
        dlgtitle = 'Dynamic data';
        dims = [1 50];
        definput = {'1000 1160 1250 1350 1550 1700 2200','300 310 325 340 367 380 390','5 5 5 5 5 5 5','100 100 100 100 100 100 100'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        m_zS=sscanf(answer{1, 1},'%f')';         % depths
        ageS=(sscanf(answer{2, 1},'%f')*1000)';  % ages (years)
        BDin=sscanf(answer{3, 1},'%f')';         % Bioturbation depths (cm)
        abuIFA=sscanf(answer{4, 1},'%f')';       % number of IFA specie individuals per cm
   
        for i=1:1:length(m_zS)
            if  mod(m_zS(i),1)>0
                error(' Error: depths values should be integers')
            end
        end

        if range(abuIFA) ~= 0 && abu==0 % if abuIFA is not constant
            error('Error\n%s', 'The variable "abu" should be 1 or 2')
        end
    else % abu=2  Temporal abundance data
        
        % load SAR and SMLD data
        prompt = {'Enter core depths (in cm !!):','Enter Ages (kaBP):','Enter Bioturbation depth (cm):'};
        dlgtitle = 'Dynamic data';
        dims = [1 50];
        definput = {'1000 1160 1250 1350 1550 1700 2200','300 310 325 340 367 380 390','5 5 5 5 5 5 5'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        m_zS=sscanf(answer{1, 1},'%f')';         % depths (cm)
        ageS=(sscanf(answer{2, 1},'%f')*1000)';  % ages (years)
        BDin=sscanf(answer{3, 1},'%f')';         % Bioturbation depths

        % load abundance data
        prompt = {'Enter Ages (kaBP):','Enter abundance:'};
        dlgtitle = 'Dynamic SAR and SMLD';
        dims = [1 50];
        definput = {'300 310 325 340 367 380 390','5 5 5 5 5 5 5'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        ageAbu=(sscanf(answer{1, 1},'%f')*1000)'; % age abundance data (years)
        abuIFA=sscanf(answer{2, 1},'%f')';        % indv abundance
    end
     
    if abu~=0
         prompt = {'A value to scale abundance e.g., 100. As the scale parameter increases, smaller temporal changes in abundance can be modeled. However, it will be more time expensive:'};
         dlgtitle = 'Scale individuals abundance';
         dims = [1 50];
         definput = {'10'};
         answer = inputdlg(prompt,dlgtitle,dims,definput);
         scale=sscanf(answer{1, 1},'%f'); % scale parameter
    else
        scale=0;
    end

     %% 1.2b. Prepare data for bioturbation simulation
     % check inputs

     l=length(ageS);
     in= length(m_zS) == length(ageS) && length(ageS) == length(BDin) && length(BDin)== length(abuIFA);
     if in ~=1
         error('Error: The lengths of Age,Sed_r,BD and abudance variables do not match, Please insert the same number of data points for each variable')
     end

     if ageS(end)<sstart_y % add cte sed_r, abu and BD 
            m_zSe=(m_zS(length(m_zS))-m_zS(length(m_zS)-1))*((sstart_y)-ageS(length(ageS)))...
            /(ageS(length(ageS))-ageS(length(ageS)-1));
            m_zS=[m_zS,m_zS(length(m_zS)) + m_zSe];
            ageS=[ageS,sstart_y];
            abuIFA=[abuIFA,abuIFA(length(abuIFA))];
            BDin=[BDin,BDin(length(BDin))]; 
            if exist("ageAbu",'var')
                  ageAbu=[ageAbu,sstart_y];
            end
     end

     if ageS(1)>send_y  % add cte SAR,abu and BD
             m_zSe=-1*(m_zS(2)-m_zS(1))*(ageS(1)-send_y)...
            /(ageS(2)-ageS(1));
             m_zS=[m_zS(1)+ m_zSe,m_zS];
             ageS=[send_y,ageS];
             abuIFA=[abuIFA(1),abuIFA];
             BDin=[BDin(1),BDin];
              if exist("ageAbu",'var')
                  ageAbu=[send_y,ageAbu];
              end
     end

     % get depths idx per month
     cum_month=cumsum(repmat(1/12,1,length(age))) + send_y; 
     depth_ind=(interp1(ageS,m_zS, cum_month,'linear','extrap'));
     depth_ind(depth_ind <0)=depth_ind(depth_ind <0) -1; % correct -0cm vs +0cm assign diff depths (if necessary)
     depths_ind=flip(depth_ind);
     pos_depth=fix(depths_ind);  % assign same depth integerer to each discrete depth
     z=max(pos_depth):-1:min(pos_depth);  % final core depths 
          
     BD=interp1(m_zS,BDin,z,'linear','extrap'); % downcore BD (cm)
     
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
             % december (Youngest age)
             seas_indx=1:1:length(age);
             ll_res=reshape(seas_indx,12,[]);
             ll_seas=reshape(ll_res(seas,:),[],1); %index seas months
             rep_ind=rep_ind(ll_seas);
            
             % take seas data and repeat it according to abu

             depths_ind=repelem(depths_ind(ll_seas),rep_ind);
             age=repelem(age(ll_seas),rep_ind);
             pos_depth=repelem(pos_depth(ll_seas),rep_ind);
             
         else  %yearly specimens
             depths_ind=repelem(depths_ind,rep_ind);
             age=repelem(age,rep_ind);
             pos_depth=repelem(pos_depth,rep_ind);
             
         end

     else  % no dynamical abundance

         if length(seas)<12     %if seasonal preference       
             % Find indx for seas months (before bio). Note that simulation ends at
             % december 
             seas_indx=1:length(age);
             ll_res=reshape(seas_indx,12,[]);
             ll_seas=reshape(ll_res(seas,:),[],1); %index seas months
             depths_ind=depths_ind(ll_seas);
             age=age(ll_seas);
             pos_depth=pos_depth(ll_seas);
             
         end
     end
     
     %% 1.3b. Bioturbation simulation

     z=flip(unique(pos_depth));  % final core depths  
     
     if any(z<0) % if needed
         z=[z,z(end)-1];
         BD=[BD, BD(end)];
     end
     
     age_bio=repmat(age',1,run); 
     
     for r= 1:run
         for i = 1:1:length(z)
             ind= find(depths_ind >= z(i) & depths_ind <z(i)+BD(i)); % find indx of depths values inside
             rndind = ind(randperm(length(ind))); % random displacement
             age_bio(ind,r)=age_bio(rndind,r);
         end
     end
     
     %% 1.4b Get IFA time slice data

     % take IFA time slice data

     age=flip(age);
     age_bio=flip(age_bio);
     depths_ind=flip(depths_ind);
     z=flip(z);
     pos_depth=flip(pos_depth);
    
     if any(z<0) % if needed
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

     age_bio_Past=cell(length(IFApast_mage),run); % bioturbated ages downcore IFA data depths

     for ii=1:length(IFApast_mage)
         [~,idx]=min(abs(med_bio-IFApast_mage(ii))); 
         for i=1:run
             age_bio_Past{ii,i} = age_bio(s_idx(idx(i)):e_idx(idx(i)),i); 
         end
     end
        
     %% 1.5b Probability of pick n older individuals on IFA time slice
    
     clearvars -except num_p run age_bio_Past IFApast_mage 

     % Prepare resampled index data
    
     mc=1000; % repeat "picking" mc times
     months_kyrs=12*1000;
     num_bio=1:20;   % number of bioturbated individuals 
     age_offset=10;  % age offset (kyr)
     num_bio_ex=reshape(repelem(num_bio,mc,age_offset),mc,age_offset,[]);
     ageoffset_ex=reshape(repelem(1:age_offset,num_p,mc),num_p,mc,[]);
     indvp=zeros(mc,age_offset,length(num_bio),run); % "individuals picked"
       
     for r=1:run % IFA data picked 1 depth
         if length(IFApast_mage) ==1 
             mage_past = mean(cell2mat(age_bio_Past(r)));
             cage_bio_Past = cell2mat(age_bio_Past(r));
         else % more than 1 depth
             mage_past=mean(reshape(vertcat(age_bio_Past{1:length(IFApast_mage),r}),[],1));
             cage_bio_Past = reshape(vertcat(age_bio_Past{1:length(IFApast_mage),r}),[],1);
         end
        
         idx = randi(length(cage_bio_Past),num_p,mc);
         older_ind= cage_bio_Past(idx)> (ageoffset_ex*months_kyrs) + mage_past;
         n_older_ind= squeeze(sum(older_ind)) > num_bio_ex;
         indvp(:,:,:,r)=n_older_ind;
     end
     pp=100*squeeze(sum(sum(indvp,1),4)/(mc*run)); % probability matrix
end

%% 2. Final Plot

% Heatmat chart, which displays the probability of picking in a given discrete
% sample 1 to 20 individuals that precipitated their shells during time 
% slices 1-to-10 kyr older than the IFA sample time slice

figure
h=heatmap(1:20,flip(1:10),flip(pp)) ;
h.Colormap = hot;
h.CellLabelFormat = '%.3f';
h.ColorLimits = [0 100];
xlabel('nº individuals')
ylabel('Age difference (kyr)')
title('Picking Probability (%)')


clearvars -except pp

