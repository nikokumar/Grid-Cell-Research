function [ NEW_STD, NEW_STL ] = Dilute_Spikes_neuron_wise( Rat,speed_threshold,tr,cluster )

% Dilute either dark or light spike trains according.
% Doing so for each neuron seperatly and according to speed threhold segments.
rng(tr);
%% Speed thresholding
ori_dir=pwd;

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


load('data.mat'); % where speeds are
speed_l=data.task(2).tracking.MY_speed; % speeds for light experiment
speed_d=data.task(1).tracking.MY_speed; % speeds for dark experiment

speed_idc_l_A=find(speed_l>=speed_threshold); % indices of speed above thresohld, light
speed_idc_l_B=find(speed_l<speed_threshold); % indices of speed below thresohld, light
speed_idc_d_A=find(speed_d>=speed_threshold); % indices of speed above thresohld, dark
speed_idc_d_B=find(speed_d<speed_threshold); % indices of speed below thresohld, dark
%% load spike trains
load('ST_D.mat'); % dark spike trains
STD=ST; % locate it 
NEW_STD=zeros(size(STD)); % will contain the modified dark spike train


load('ST_L.mat'); % light spike trains
STL=ST;
NEW_STL=zeros(size(STL)); % will contain the modified light spike train

N=size(ST,1); % # of neurons
clear ST;

cd(ori_dir); % going back to directory we came from


STL_Above = STL(:,speed_idc_l_A);
STL_Below = STL(:,speed_idc_l_B);
STD_Above = STD(:,speed_idc_d_A);
STD_Below = STD(:,speed_idc_d_B);
%% Compare neurons wise
% ld is a marker lableing which spiek train need dilution:
% ld=0: dilute light spikes, ld=1: dilute dark spikes

for k=1:N % running on all neurons
    
    sl_A=STL_Above(k,:); % spike train in times above speed threshold for light
    sl_B=STL_Below(k,:); % spike train in times below speed threshold for light
    
    sd_A=STD_Above(k,:); % spike train in times above speed threshold for dark
    sd_B=STD_Below(k,:); % spike train in times below speed threshold for dark
    
    for S=1:2 % running on above and below segments
        if S==1 % analyze below segment
            SD=sd_B; % dark below
            SL=sl_B; % light below
        end
        if S==2 % analyze above segment
            SD=sd_A; % dark above
            SL=sl_A; % light above
        end
    %% Diluting above speed threshold segment
        pf = mean(SD)/mean(SL); % ratio between dark and lights mean f.r for times above threshold
        st=SL; % the spike train to dilute, we dilute the light spike trains by defult
        ld=0; % mark that we dilute light spike train, by defult
        if pf>=1 % if the are more mean spikes in dark (for times above speed threhold)
            st=SD; % the spike train to dilute, actually we need to dilute the dark spike train
            ld=1; % mark that we actually dilute dark spike train
            pf=1/pf; % the percentage to use
        end

        per = 1-pf; % percent to use for the dilution


        Maximum_spikes=max(st); % maximal numer of spikes per time bin
        spikes_locs=cell(1,Maximum_spikes);
        for j=1:Maximum_spikes
            spikes_locs{j}=find(st>j-1)';
        end

        Spikes_loc=[]; % will contain all spikes locations
        for j=1:Maximum_spikes
            Spikes_loc=[Spikes_loc;spikes_locs{j}]; % unrolling all spikes 
        end


        fil=rand(size(Spikes_loc));
        fil(fil<=per)=1; % spikes to filter
        fil(fil~=1)=0;

        reduce_spike=Spikes_loc(fil==1); % spikes to filter
        if numel(unique(reduce_spike))==numel(reduce_spike) % if all values are unique we can do it in a ingle line
            st(reduce_spike)=st(reduce_spike)-1; % diluting
        end

        if numel(unique(reduce_spike))~=numel(reduce_spike) % if not than we nee to do a loop...
            for j=1:numel(reduce_spike) % we do it in a loop and not in 1 vector iteration because of non unique values
                st(reduce_spike(j))=st(reduce_spike(j))-1; % diluting
            end
        end
    
        %% Allocating and building back modified spike trains
        if S==1 % below case
            if ld==0 % if we diluted the light spike train we modify it, and dark stays unchanged
                
                NEW_STD(k,speed_idc_d_B) = STD(k,speed_idc_d_B); % dark unchanged, same as original in the below segments
                NEW_STL(k,speed_idc_l_B) = st; % modified light spike train, in the below segments

            elseif ld==1 % if we diluted the dark spike train we modify it, and light stays unchanged
                
                NEW_STD(k,speed_idc_d_B) = st; % modified dark spike train, in the below segments
                NEW_STL(k,speed_idc_l_B) = STL(k,speed_idc_l_B); % light unchanged, same as original in the below segments
                
            end
        end
        
        if S==2 % above case
            if ld==0 % if we diluted the light spike train we modify it, and dark stays unchanged
                    
                NEW_STD(k,speed_idc_d_A) = STD(k,speed_idc_d_A); % dark unchanged, same as original in the above segments
                NEW_STL(k,speed_idc_l_A) = st; % modified light spike train, in the above segments
                
            elseif ld==1 % if we diluted the dark spike train we modify it, and light stays unchanged
                
                NEW_STD(k,speed_idc_d_A) = st; % modified dark spike train, in the above segments
                NEW_STL(k,speed_idc_l_A) = STL(k,speed_idc_l_A); % light unchanged, same as original in the above segments
            end
        end
    end
end

end