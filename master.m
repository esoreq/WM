%% Working memory results script for submission
%% Data prepration

% create.Data; run once to copy all data 
%% Behavioural results
% Study 1/2 accuracy & Study 1/2 RT
%%
create.Behavioural; % 
%% Mapping multiple demand areas
%  Conjunction results (DG & SG maps) across studies
% 
% And cross study corrospondence eproducibility across studies
%%
create.Conjunction;
%% ROI definition with Watershed algorithm
%%
create.Parcelation;
%% Hierarchical cluster analysis using activation
%%
create.HCL;
%% Maintenance through connectivity: global analyses
%%
create.GlobalEffects;
%% Multinomial classification of factors (SVM ECOC)
%%
create.Decoding;

%% Projection of sparse models onto the brain
%%
create.ModelsMNIprojection;