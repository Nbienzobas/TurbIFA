function results=quantifaerrv3(X,Y,anerr,exerr,num_q,X_tslice,Y_tslice,run,varargin)

% About
% *************************************************************************
%
% quantifaerrv3 is a function that calculates the i) noise (IFA size + analitycal) 
% in IFA datasets 1σ value and in a quantile-quantile space and,
% ii) the probability that a downcore change (increase/decrease) in the IFA-1σ 
% is significant
%
% --> IFA distribution 1-σ value: quantifaerr builds on INFAUNAl algorithm
%     to calculate the probability of a significant change (95%) in IFA distribution
%     1-σ value. It uses the root-mean square error of the IFA core-top and 
%     downcore standard errors as a threshold for significance
%
% --> IFA distribution in a q-q space: quantifaerr builds on QUANTIFA algorithm  
%     to calculate the uncertainties (95% confidence interval) associated with 
%     contrasting the core-top vs downcore IFA distributions in the q-q space 
%
% Source INFAUNAL
% https://www.mathworks.com/matlabcentral/fileexchange/47695-infaunal
% https://doi.org/10.1002/palo.20037
%
% Source QUANTIFA: 
% https://github.com/rh-glaubke/QUANTIFA
% https://doi.org/10.5281/zenodo.7775163
% https://doi.org/10.1029/2020PA004065
%
%
% Syntax: quantifaerrv3(X,Y,anerr,num_q,X_tslice,Y_tslice,run,varargin)
%
% Inputs:
% *************************************************************************
% X = reference IFA dataset
% Y = other IFA dataset
% anerr = anaytical uncertainty
% exerr = extra source of gaussian noise
% num_q = number of quantiles
% mc = number of times picking process is repeated
% X_tslice = X time series to resampled
% Y_tslice = Y time series to resampled
% run = number of bioturbated synthethic cores
%
% Main Outputs:
% *************************************************************************
% results.Xerrqq = uncertainty X IFA dataset (as percentiles 2.5 & 97.5, q-q space)  
% results.Yerrqq = uncertainty Y IFA dataset (as percentiles 2.5 & 97.5, q-q space)
% results.Xstd_per = uncertainty X IFA dataset 1-sigma value (as percentiles 2.5 & 97.5)
% results.Ystd_per = uncertainty Y IFA dataset 1-sigma value (as percentiles 2.5 & 97.5)
% results.std_i  = probability of a significant downcore increase in IFA
%                  distribution 1-sigma value
% results.std_d  = probability of a significant downcore decrease in IFA
%                  distribution 1-sigma value
%
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 5-Jan-2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Declare data

ng=nargin;
if ng==9
    fig_24=varargin{1};  % figures 2-4 ms
elseif ng==10             % bioturbated data
    IFAtop_mage=varargin{1};
    IFApast_mage=varargin{2};
end
mc=1000; % pick IFA datasets 1000 times
rx=[length(X),mc,2]; % size "picked" X population dataset
ry=[length(Y),mc,2]; % size "picked" Y population dataset

if iscell(X_tslice)== 0  

     % X population "picking" + error terms   
     mc_x=X_tslice(randi(length(X_tslice),rx)) + anerr*randn(rx) + ...
         exerr*randn(rx);
     
     % Y population "picking" +  error terms     
     mc_y=Y_tslice(randi(length(Y_tslice),ry)) + anerr*randn(ry) + ...
         exerr*randn(ry);

     % Uncertainty IFA 1σ value (bioturbated populations)
     mc_xx=squeeze(mc_x(:,:,1));
     mc_yy=squeeze(mc_y(:,:,1));
     % X population
     sigma_X= std(X)+ (std(mc_xx)- mean(std(mc_xx))); 
     per_sX=prctile(std(mc_xx)- mean(std(mc_xx)),[2.5,97.5]); % 95
     % Y population
     sigma_Y= std(Y) + (std(mc_yy)-mean(std(mc_yy))); 
     per_sY=prctile(std(mc_yy)- mean(std(mc_yy)),[2.5,97.5]); % 95 
     % calculate diff 
     difsig= sigma_Y - sigma_X';
     % prob significant downcore increase IFA 1-std value
     incres= difsig>0 & difsig > sqrt(per_sY(1)^2 + per_sX(2)^2);
     sig_i= 100*(sum(sum(incres))/(mc*mc));
     % prob significant downcore decrease IFA 1-std value
     decreas= difsig<0 & abs(difsig) > sqrt(per_sY(2)^2 + per_sX(1)^2);
     sig_d= 100*(sum(sum(decreas))/(mc*mc));
     sigmaX_per=per_sX;
     sigmaY_per=per_sY;
    
else 
   
    % X population "picking" + error terms                          
    mc_x= zeros(length(X),mc,run,2); 
     
    for i=1:1:run
        Xbio_i=reshape(vertcat(X_tslice{1:length(IFAtop_mage),i}),[],1);
        mc_x(:,:,i,:) = Xbio_i(randi(length(Xbio_i),rx)) + ...
            anerr*randn(rx) + exerr*randn(rx); 
    end 

    % Y population "picking" + error terms 
    mc_y= zeros(length(Y),mc,run,2); 
  
    for i=1:1:run
        Ybio_i=reshape(vertcat(Y_tslice{1:length(IFApast_mage),i}),[],1);
        mc_y(:,:,i,:) = Ybio_i(randi(length(Ybio_i),ry)) + ...
            anerr*randn(ry) + exerr*randn(ry); 
    end

    % Uncertainty IFA 1σ value (bioturbated populations)
    sig_i=zeros(run,1);  % prob significant downcore increase IFA 1-std value
    sig_d=zeros(run,1);  % prob significant downcore decrease IFA 1-std value
    sigmaX_per=zeros(2,run);
    sigmaY_per=zeros(2,run);

    for i=1:run
        mc_xx=squeeze(mc_x(:,:,i,1));
        mc_yy=squeeze(mc_y(:,:,i,1));
        % X population
        sigma_X= std(X)+ (std(mc_xx)- mean(std(mc_xx))); 
        per_sX=prctile(std(mc_xx)- mean(std(mc_xx)),[2.5,97.5]); % 95
        % Y population
        sigma_Y= std(Y) + (std(mc_yy)-mean(std(mc_yy))); 
        per_sY=prctile(std(mc_yy)- mean(std(mc_yy)),[2.5,97.5]); % 95 
        % calculate diff 
        difsig= sigma_Y - sigma_X';
        % prob significant downcore increase IFA 1-std value
        incres= difsig>0 & difsig > sqrt(per_sY(1)^2 + per_sX(2)^2);
        sig_i(i)= 100*(sum(sum(incres))/(mc*mc));
        % prob significant downcore decrease IFA 1-std value
        decreas= difsig<0 & abs(difsig) > sqrt(per_sY(2)^2 + per_sX(1)^2);
        sig_d(i)= 100*(sum(sum(decreas))/(mc*mc));
        sigmaX_per(:,i)=per_sX;
        sigmaY_per(:,i)=per_sY;
    end
        
end


if ndims(mc_x)<4 % non bioturbated data or results for figures 2-4 ms
    
    % tranform data to q-q space
    if  exist('fig_24','var')==1 % figures 2_4 ms (non-normalized qq)
            mc_qx = quantile(mc_x,num_q); % X population
            mc_qy= quantile(mc_y,num_q);  % Y population
    else
           mc_qx = quantile(mc_x-mean(mc_x),num_q); % X population
           mc_qy = quantile(mc_y-mean(mc_y),num_q); % Y population
    end
    

    % ---- Generate Bins ----

    % X population

    bx = mean(mc_qx(:,:,1),2);  
    binsx = quantile(bx,num_q)';                
    dx = diff(binsx)/2; 

    % Y population

    by = mean(mc_qy(:,:,1),2);  
    binsy = quantile(by,num_q)';                
    dy = diff(binsy)/2; 

    edges = zeros(length(binsx)+1,2);  
   
   for ii = 1:2
       for i = 2:length(binsx)
           edges(i,1) = binsx(i-1) + dx(i-1);
           edges(i,2) = binsy(i-1) + dy(i-1);
       end
       if ii == 1
           edges(1,ii) = min(mc_qx(1,:,1));
           edges(end,ii) = max(mc_qx(num_q,:,1));
       else
           edges(1,ii) = min(mc_qy(1,:,1));
           edges(end,ii) = max(mc_qy(num_q,:,1));
       end
   end
   
   % ---- Bin Data ----
   serrX=zeros(2,num_q); % uncertainty for X population q-q space
   serrY=zeros(2,num_q); % uncertainty for Y population q-q space
   
   for ii = 1:2                        % Calculate uncertainy for X data
       if ii == 1
           for i = 1:num_q
               if i < num_q
                   idx = find(mc_qx(:,:,1)>=edges(i,ii) ...
                    & mc_qx(:,:,1)<edges(i+1,ii));
               else
                   idx = find(mc_qx(:,:,1)>=edges(i,ii) ...
                    & mc_qx(:,:,1)<=edges(i+1,ii));
               end
               A = mc_qx(:,:,2);
               qx2_vals = A(idx);
               per = prctile(qx2_vals-mean(qx2_vals), [2.5, 97.5]);
               serrX(:,i) = per';
           end
       else
           for i = 1:num_q             % Calculate uncertainty for Y data
               if i < num_q
                   idx = find(mc_qy(:,:,1)>=edges(i,ii) ...
                    & mc_qy(:,:,1)<edges(i+1,ii));
               else
                   idx = find(mc_qy(:,:,1)>=edges(i,ii) ...
                    & mc_qy(:,:,1)<=edges(i+1,ii));
               end
               A = mc_qy(:,:,2);
               qy2_vals = A(idx);
               per = prctile(qy2_vals-mean(qy2_vals), [2.5, 97.5]);
               serrY(:,i) = per';
             
           end
       end
   end
 
else  % bioturbated datasets + n runs mixing model

     % Declare data
     serrX=zeros(2,num_q,run); % uncertainty for X population q-q space
     serrY=zeros(2,num_q,run); % uncertainty for Y population q-q space
     qxr=zeros(num_q,mc,2,run);
     qyr=zeros(num_q,mc,2,run);

    for r=1:run
        mc_xr= squeeze(mc_x(:,:,r,:));
        mc_yr= squeeze(mc_y(:,:,r,:));

        mc_qx=zeros([num_q,size(mc_xr,2),2]);
        mc_qy=zeros([num_q,size(mc_yr,2),2]);
       
       for i=1:2
           mc_qx(:,:,i) = quantile(mc_xr(:,:,i)-mean(mc_xr(:,:,i)),num_q); % X population
           mc_qy(:,:,i) = quantile(mc_yr(:,:,i)-mean(mc_yr(:,:,i)),num_q); % Y population
       end
       
       % ---- Generate Bins ----

       % X population

       bx = mean(mc_qx(:,:,1),2);  
       binsx = quantile(bx,num_q)';                
       dx = diff(binsx)/2; 

       % Y population

       by = mean(mc_qy(:,:,1),2);  
       binsy = quantile(by,num_q)';                
       dy = diff(binsy)/2; 

       edges = zeros(length(binsx)+1,2);

       for ii = 1:2
           for i = 2:length(binsx)
               edges(i,1) = binsx(i-1) + dx(i-1);
               edges(i,2) = binsy(i-1) + dy(i-1);
           end
           if ii == 1
               edges(1,ii) = min(mc_qx(1,:,1));
               edges(end,ii) = max(mc_qx(num_q,:,1));
           else
               edges(1,ii) = min(mc_qy(1,:,1));
               edges(end,ii) = max(mc_qy(num_q,:,1));
           end
       end
       
       % ---- Bin Data ----
      
       
       for ii = 1:2  % Calculate uncertainy for X data                      
           if ii == 1
               for i = 1:num_q
                   if i < num_q
                       idx = find(mc_qx(:,:,1)>=edges(i,ii) ...
                       & mc_qx(:,:,1)<edges(i+1,ii));
                   else
                       idx = find(mc_qx(:,:,1)>=edges(i,ii) ...
                       & mc_qx(:,:,1)<=edges(i+1,ii));
                   end
                   A = mc_qx(:,:,2);
                   qx2_vals = A(idx);
                   per = prctile(qx2_vals-mean(qx2_vals), [2.5, 97.5]);
                   serrX(:,i,r) = per';
               end
           else
               for i = 1:num_q    % Calculate uncertainty for Y data         
                   if i < num_q
                       idx = find(mc_qy(:,:,1)>=edges(i,ii) ...
                       & mc_qy(:,:,1)<edges(i+1,ii));
                   else
                       idx = find(mc_qy(:,:,1)>=edges(i,ii) ...
                       & mc_qy(:,:,1)<=edges(i+1,ii));
                   end
                   A = mc_qy(:,:,2);
                   qy2_vals = A(idx);
                   per = prctile(qy2_vals-mean(qy2_vals), [2.5, 97.5]);
                   serrY(:,i,r) = per';
               end
           end
       end
       qxr(:,:,:,r)=mc_qx;
       qyr(:,:,:,r)=mc_qy;
    end
end

% save results

if iscell(X_tslice)== 0   
   results.mc_x=mc_x;   % resampled Y population data
   results.mc_y=mc_y;   % resampled Y population data
end

if exist('sigmaX_per','var')==1  
    results.Xstd_per=sigmaX_per; % uncertainty X population spread (1σ distribution value)
    results.Ystd_per=sigmaY_per; % uncertainty Y population spread (1σ distribution value)
    results.std_i=sig_i;
    results.std_d=sig_d;
end

if ndims(mc_x)<4
    results.qqx=mc_qx(:,:,1);
    results.qqy=mc_qy(:,:,1);
else
    results.qqx=qxr(:,:,1,:);
    results.qqy=qyr(:,:,1,:);
end

results.Xqq_err=serrX;  % uncertainty X population qq space
results.Yqq_err=serrY;  % uncertainty Y population qq space
