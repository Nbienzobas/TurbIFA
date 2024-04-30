function ageBH=analiBH(sed_r,BD,nIFA)
%
% function to calculate the probability of pick 1-20 individuals 
% 1-10 kyrs older than (mean) age of a cm
% Solves analitically Berger and Heath model by calculating and sampling 
% an exponential age distribution with an scale parameter equal to SMLD/SAR
%
% Author: Natalia Bienzobas Montávez 
% Centro de Investigación Mariña, Universidade de Vigo, GEOMA,
% Palaeoclimatology Lab, Vigo, 36310,Spain
% email addresses: nbienzovas@uvigo.gal
% Last revision: 20-Dec-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_bio=1:20;   % number of bioturbated individuals 
age_offset=10;  % age offset (kyr)
mc=1000;        % "pick" 1000 times
kyr=1000;
sizedis=10000000; % to force the distribution to fit B&H solution
num_bio_ex=reshape(repelem(num_bio,mc,age_offset),mc,age_offset,[]);
ageoffset_ex=reshape(repelem(1:age_offset,length(nIFA),mc),length(nIFA),mc,[]);
scale=BD(1)/sed_r;
% exponential random distribution
age_dis=(-scale*log(rand(1,sizedis)))*kyr; 
% picking
idx = randi(length(age_dis),length(nIFA),mc);  
ind_sam= age_dis(idx)> (ageoffset_ex*kyr) + mean(age_dis);
oldf= squeeze(sum(ind_sam)) > num_bio_ex;
% probability
pbt=100*squeeze(sum(sum(oldf,1),4)/(mc)); 
ageBH.pp=pbt;
end