function [ST_Per ] = permuted_spike_train( Rat,ld,ST,modules_together,speed_threshold,tr)
% permute (reshuffle) the spiking locations for each neuorn to generate a
% new spike train and calculate its likelihood - compare to the likelihood
% from pure dark trials without manupulation

% inputs: (1) ST: = the original spike train
%         (2) ld: light or dark (0 for dark and 1 for light)
%         (3) modules_together: if ==1 then we apply the same
%         permutation to all neurons that belong to the same modle, if not
%         than we apply random independent permutation to each neuron
%         (4) speed_threshold: speed threshold
%         
%  output: ST_Per = permuted spike train

N=size(ST,1); % # of neurons
T=size(ST,2); % # of time bins
ST_Per=zeros(N,T); % will contain the final new permuted spike trains for all neurons

ori_dir=pwd;
if Rat==1 % Bubble
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble/Rat_Bubble_data'); % cluster run
%     cd('Rat_Bubble_data'); % local run
elseif Rat==11 % Bubble 2nd dataset
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Bubble 2/Rat_Bubble_data 2'); % cluster run
%     cd('Rat_Bubble_data 2'); % local run
elseif Rat==2 % Roger
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Roger/Rat_Roger_data'); % cluster run
%     cd('Rat_Roger_data'); % local run
elseif Rat==3 % Coconut
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Coconut/Rat_Coconut_data'); % cluster run
%     cd('Rat_Coconut_data'); % local run
elseif Rat==4 % Rolf
    cd('/ems/elsc-labs/burak-y/haggai.agmon/Torgeir New data Rolf/Rat_Rolf_data'); % cluster run
%     cd('Rat_Rolf_data'); % local run
end

rat_dir=pwd;
load('data.mat'); % where speeds are
cd(ori_dir);
%%
speed=data.task(ld+1).tracking.MY_speed; % speeds
speed_idc_A=find(speed>=speed_threshold); % indices of speed above thresohld
speed_idc_B=find(speed<speed_threshold); % indices of speed below thresohld

ST_A = ST(:,speed_idc_A); % spike trains for times above threshold
ST_B = ST(:,speed_idc_B); % spike trains for times below threshold
T_A = size(ST_A,2); % # of time bins above threshold
T_B = size(ST_B,2); % # of time bins below threshold

tr=(tr-1)*(N+1); % random seed generator
%% Same permutation to all neurons
rng(tr);
% Permutation_A=randperm(T_A); % the permutation that will be applied to the k'th neuron
% Permutation_B=randperm(T_B); % the permutation that will be applied to the k'th neuron
%% Independent permutations to all neurons
if modules_together~=1 % meaning each neuron permutation is performed independetly

    for k=1:N % running on each neurons sepratly because there is not enough memory to do them all together

        st_a=ST_A(k,:); % original spike train of the k'th neuron on segments above thresohld
        st_b=ST_B(k,:); % original spike train of the k'th neuron on segments below thresohld

        %% Different permutation to each neuron
        rng(k+tr); % can make the same permutation (or a different one) to all the neurons
        %% Continue
        Permutation_A=randperm(T_A); % the permutation that will be applied to the k'th neuron for above threshold segments
        Permutation_B=randperm(T_B); % the permutation that will be applied to the k'th neuron for below threshold segments
       
        ST_Per(k,speed_idc_A) = st_a(Permutation_A); % the new permuted spike train of the k'th neuron for above threshold segments
        ST_Per(k,speed_idc_B) = st_b(Permutation_B); % the new permuted spike train of the k'th neuron for below threshold segments
    end   
end

%% Same permutation to all neurons that belong to the same module
if modules_together==1
    cd(rat_dir);
    load('Neurons_sorted_according_to_modules.mat');
    load('ID.mat');
    cd(ori_dir);
    
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
    
    M_ID=cell(1,3);
    M_ID{1}=m1_id;
    M_ID{2}=m2_id;
    M_ID{3}=m3_id;

    for m=1:numel(M_ID) % running on all modules
        m_id=M_ID{m}; % all neurons that belong to this module
        
        rng(m+tr); % can make the same permutation (or a different one) to all the neurons
      
        Permutation_A=randperm(T_A); % the permutation that will be applied to the k'th neuron for above threshold segments
        Permutation_B=randperm(T_B); % the permutation that will be applied to the k'th neuron for below threshold segments
        
        for k=1:numel(m_id) % running on all neurons that belong to the m'th modules 
            
            h=m_id(k); % the relevant neuron from this module
            
            st_a=ST_A(h,:); % original spike train of the k'th neuron on segments above thresohld
            st_b=ST_B(h,:); % original spike train of the k'th neuron on segments below thresohld
     
            ST_Per(h,speed_idc_A) = st_a(Permutation_A); % the new permuted spike train of the k'th neuron for above threshold segments
            ST_Per(h,speed_idc_B) = st_b(Permutation_B); % the new permuted spike train of the k'th neuron for below threshold segments
        end
    end
end