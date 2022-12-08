function [  ] = Pair_wise_correlations_for_Synthetic_spikes(Rat,L_D,L,Groups,G,tr )
% Correlations between ineter and intra module pair-wise correlations
% 1. Generate synthetic data based on recorded data
% 2. Measure it's cross correlation as a function of the arena size


cluster=0; % run on cluster? set to 1

% Rat=1; % 1 for Bubble, 2 for Roger
% L_D=0; % analyze dark or light data? 0 for dark and 1 for light
% speed_threshold=0; % speed threshold

max_lag=2; % maximum time lag [sec]

if Groups==1 % if we plot scatter plot, we look only at the zero lag correlations
    max_lag=0;
end

if L_D==0
    Na='Dark';
elseif L_D==1
    Na='Light';
end

ori_dir=pwd; % current directory
%% Check if result already exist on cluster
if cluster==1
    
    cd(Na);
    r_dir=pwd;
    
    name = [Na,sprintf(':Radi=%dcm,Groups=%d,G=%d,tr=%d.mat',L,Groups,G,tr)]; % name of file
    if exist(name,'file')~=0 % if this file exists we don't run the script
        return;
    end

    cd(ori_dir); % going back to the directory we came from
end
%% Loading data
if Rat==1 % Bubble
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Rat_Bubble_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data'); % when run locally
    end
elseif Rat==11 % Bubble 2nd dataset
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2/Rat_Bubble_data 2'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data 2'); % when run locally
    end
elseif Rat==2 % Roger
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger/Rat_Roger_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Roger_data'); % when run locally
    end
elseif Rat==3 % Coconut
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut/Rat_Coconut_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Coconut_data'); % when run locally
    end
elseif Rat==4 % Rolf
    if cluster==1
        cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf/Rat_Rolf_data'); % when run on cluster
    elseif cluster~=1
        cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Rolf_data'); % when run locally
    end
end
% rat_dir=pwd;

load('data.mat'); % loading the data file
Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]

lag_ind = round(max_lag/dt); % how many indices are needed for the time lag

load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
FN=fieldnames(Neurons_sorted_according_to_modules); % names of fields
Modules = numel(FN); % # of modules

N=zeros(1,Modules); % will contain the # of neurons in each module

for nn=1:Modules
    N(nn) = numel(Neurons_sorted_according_to_modules.(FN{nn})); % # of neurons in the nn'th module
end
cd(ori_dir); % go back to original directory where we came from
%% Allocation to modules
% Corresponding allocation to modules of the neurons which their spikes are
% in ST file and tuning curve 
m1_id = 1:N(1);
m2_id = m1_id(end)+1:m1_id(end)+N(2);
m3_id = m1_id(end)+N(2)+1:sum(N);
%% Generate synthetic spike trains
[STs,TCs] = Sythetic_spikes(Rat,Modules,N,L_D,dt,L,tr,cluster);

ST = [STs{1} ; STs{2}; STs{3}]; % spike trains of all neurons
clear STs;

Tuning_Curves = zeros(2*L+1,2*L+1,sum(N));
Tuning_Curves(:,:,m1_id)=TCs{1};
Tuning_Curves(:,:,m2_id)=TCs{2};
Tuning_Curves(:,:,m3_id)=TCs{3};
clear TCs;
%% Grouping neurons
g1=numel(m1_id)/Groups; % # of neurons in each group for 1st module
g2=numel(m2_id)/Groups; % # of neurons in each group for 2nd module
g3=numel(m3_id)/Groups; % # of neurons in each group for 3rd module

m1_id = m1_id(round(G*g1) +1 : round((G+1)*g1)); % neurons from 1st module in the G group
m2_id = m2_id(round(G*g2) +1 : round((G+1)*g2)); % neurons from 2nd module in the G group
m3_id = m3_id(round(G*g3) +1 : round((G+1)*g3)); % neurons from 3rd module in the G group
%% Smooth data using a gaussian kernel
sigma = .05; % desired sigma of gaussian in sec [50msec]
sigma = round(sigma/dt); % sigma in indices (dt units)

fac=5; % factor between std and window size of smoothdata gaussian function
win = fac*sigma; % the window given as an input to the smoothdata function

[smoothed_ST,~] = smoothdata(ST,2,'gaussian',win);
clear ST;
%% Intra module pairs
Intra_Pairs = cell(1,Modules); % will contain the intra module pairs

Intra_Pairs{1} = nchoosek(m1_id,2); % pairs within module 1
Intra_Pairs{2} = nchoosek(m2_id,2); % pairs within module 2
if Rat~=3 % if it's not Coconut we have a third module
    Intra_Pairs{3} = nchoosek(m3_id,2); % pairs within module 3
    Intra_Pairs = [Intra_Pairs{1}; Intra_Pairs{2}; Intra_Pairs{3}];
elseif Rat==3 % if it's Coconut we have only 2 modules
    Intra_Pairs = [Intra_Pairs{1}; Intra_Pairs{2}];
end

n_intra_pairs = size(Intra_Pairs,1); % # of intra pairs
%% Inter module pairs
Inter_Pairs = cell(1,factorial(Modules)/2); % will contain the inter module pairs

[m,n]=ndgrid(m1_id,m2_id);
Inter_Pairs{1} = [m(:),n(:)]; % pairs between modules 1 and 2

if Rat~=3 % if it's not Coconut we have 2 more pairs
    [m,n]=ndgrid(m1_id,m3_id);
    Inter_Pairs{2} = [m(:),n(:)]; % pairs between modules 1 and 3

    [m,n]=ndgrid(m2_id,m3_id);
    Inter_Pairs{3} = [m(:),n(:)]; % pairs between modules 2 and 3

    Inter_Pairs = [Inter_Pairs{1}; Inter_Pairs{2}; Inter_Pairs{3}];
elseif Rat==3 % if it's Coconut we have only 2 modules
    Inter_Pairs = [Inter_Pairs{1}];
end

n_inter_pairs = size(Inter_Pairs,1); % # of inter pairs
%% Sort the inter module pairs based on overlap between their tuning curves
overlap_inter=zeros(1,n_inter_pairs); % will contain the overlaps for each pair of neurons (from the same module)
for k=1:n_inter_pairs % running on all neurons form the same module
    overlap_inter(k)=sum(sum(Tuning_Curves(:,:,Inter_Pairs(k,1)).*Tuning_Curves(:,:,Inter_Pairs(k,2)))); % overlap of a pair
end
    
[~,I] = sort(overlap_inter); % sorting neuron in ascending order
Inter_Pairs=Inter_Pairs(I,:); % inter pairs sorted in ascending oreder according to overlap in their receptive fields
%% Cross-correlations
CORR_intra = zeros(n_intra_pairs,2*lag_ind+1); % will contain cross corr of intra pairs
MY_CORR_intra=CORR_intra; % using my normalization
CORR_inter = zeros(n_inter_pairs,2*lag_ind+1); % will contain cross corr of inter pairs
MY_CORR_inter=CORR_inter; % using my normalization

for n=1:n_intra_pairs
    MY_CORR_intra(n,:) = My_Cross_Correlation(smoothed_ST(Intra_Pairs(n,1),:) , smoothed_ST(Intra_Pairs(n,2),:),lag_ind,0);
    CORR_intra(n,:) = xcorr(smoothed_ST(Intra_Pairs(n,1),:) , smoothed_ST(Intra_Pairs(n,2),:),lag_ind,'normalized');
end
for n=1:n_inter_pairs
    MY_CORR_inter(n,:) = My_Cross_Correlation(smoothed_ST(Inter_Pairs(n,1),:) , smoothed_ST(Inter_Pairs(n,2),:),lag_ind,0);
    CORR_inter(n,:) = xcorr(smoothed_ST(Inter_Pairs(n,1),:) , smoothed_ST(Inter_Pairs(n,2),:),lag_ind,'normalized');
end
%% Organizing Results
Corrs.MY_CORR_intra=MY_CORR_intra;
Corrs.CORR_intra=CORR_intra;
Corrs.MY_CORR_inter=MY_CORR_inter;
Corrs.CORR_inter=CORR_inter;

%% Saving
cd(r_dir);
save(name,'Corrs');
end