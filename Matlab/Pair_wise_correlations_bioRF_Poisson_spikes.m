function [  ] = Pair_wise_correlations_bioRF_Poisson_spikes(Rat,L_D,speed_threshold,Groups,G,tr)
% Correlations between ineter and intra module pair-wise correlations
% Following the torus paper:
% For recorded tracjectory and tuning curves but generating Poisson spikes

rng(tr); % new realization per Poisson sampling
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
    
    name = [Na,sprintf(':speed_threshold=%d,Groups=%d,G=%d,tr=%d.mat',speed_threshold,Groups,G,tr)]; % name of file
    if exist(name,'file')~=0 % if this file exists we don't run the script
        return;
    end

    cd(ori_dir); % going back to the trajectory we came from
end
%% Loading data
Modules=3; % defult # of modules

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
    Modules=2; % total # of modules for Coconut
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

load('Tuning_Curves.mat'); % loading the tuning curves
load('ID.mat'); % corrseponding ID's of neurons
load('Neurons_sorted_according_to_modules.mat'); % loading sorting of neurons according to modules
N=size(Tuning_Curves,3); % # of total neurons (grid cells)

rounding=0; % [cm], we will work in 1cm unit resulotion

X_traj =roundn( data.task(L_D+1).tracking.x, rounding); % x location of the animal's trajectory [cm]
% rounded to 1cm resulotion
Y_traj = roundn( data.task(L_D+1).tracking.y, rounding); % y location of the animal's trajectory [cm]
% rounded to 1cm resulotion

speed=data.task(L_D+1).tracking.MY_speed; % recorded speed

lag_ind = round(max_lag/dt); % how many indices are needed for the time lag
%% Generating Poisson spikes
ST=zeros(N,numel(X_traj)); % will contain the spikes vs time
 
L=75; % [cm] radius of the circled arena
X_traj=X_traj+L+1; % set indices to correspond to rf matrix
Y_traj=Y_traj+L+1; % set indices to correspond to rf matrix

% dr=1;
% [X,Y] = meshgrid(-L:dr:L, -L:dr:L);
% R_sq=X.^2 + Y.^2; % squared distance from origin
% [out_y,out_x]=find(R_sq>L^2); % coordinates to set to 0 firing rate since they are outsile of circle
% idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

for j=1:N % running on all neurons in that module
    rf = Tuning_Curves(:,:,j); % receptive field of the j'th neuron
%     rf(idc)=0; % Sanity check - 0 firing rate outside the arena (the trajectory can't also exceed the arena boundaries)

    idx = sub2ind(size(rf),Y_traj,X_traj); % translate trajectory long indices to short index
    fr_vs_time = rf(idx); % firing rate vs time

    ST(j,:) = poissrnd(fr_vs_time.*dt); % Poisson spikes       
end
%% Allocation to modules
m1_id=Neurons_sorted_according_to_modules.M1; % id of neuorons from first module
m2_id=Neurons_sorted_according_to_modules.M2; % id of neuorons from second module
m3_id=Neurons_sorted_according_to_modules.M3; % id of neuorons from third module

for m=1:numel(m1_id) % running on 1st module neurons
    m1_id(m)=find(m1_id(m)==ID); % translating unit number to serial neuron # (from units # to 1 to N)
end
for m=1:numel(m2_id) % running on 2nd module neurons
    m2_id(m)=find(m2_id(m)==ID);
end
for m=1:numel(m3_id) % running on 3rd module neurons
    m3_id(m)=find(m3_id(m)==ID);
end
cd(ori_dir); % go back to original directory where we came from

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

if speed_threshold==1 % setting speed threshold
    speed_threshold=3;
end
smoothed_ST = smoothed_ST(:,speed>=speed_threshold); % thresholding speed keeping whats above the threshold
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