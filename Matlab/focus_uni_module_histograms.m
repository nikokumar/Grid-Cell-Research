close all; clear;
ori_dir=pwd;

Rat=11; % which rat?
tau=100; % exp time constant [ms]
L_D=0; % Dark/Light
speed_thresh=3; % speed threshold used in dilution function

thresh=.8; % include data up to covering thresh of the 2-d histogram (80%)

delta_err_multi = 3; % MAE bin width
Modules=1:3; % # of modules
if Rat==3
    Modules=1:2; % # of modules for Coconut
end
mo_pairs=combntns(Modules,2); % possible combinations of pairs
pairs=factorial(numel(Modules))/2; % # of pairs of modules

if L_D==1
    Na='Light';
elseif L_D==0
    Na='Dark';
end
%% Getting speed
if Rat==1 % Bubble
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data'); % going into rat 1 (Bubble) data directory
    ra='Bubble';
elseif Rat==11 % Bubble 2nd dataset
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data 2'); % going into rat 1 (Bubble) 2nd data directory
    ra='Bubble 2';
elseif Rat==2 % Roger
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Roger_data'); % going into rat 2 (Roger) data directory
    ra='Roger';
elseif Rat==3 % Coconut
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Coconut_data'); % going into rat 3 (Coconut) data directory
    ra='Coconut';
elseif Rat==4 % Rolf
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Rolf_data'); % going into rat 4 (Cocnut) data directory
    ra='Rolf';
end
load('data.mat'); % loading the data file
%% Multi module error to true location and mle
cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Exp decoder');
cd(ra);
cd('0 shifts results for focus');
r_dir=pwd;
cd('Multi');

load([Na,sprintf(':Tau=%dms,spatial_shift=0,tr=1.mat',tau)]);
Multi_err_to_true = Simulation_Results.ERR; % dark error vs time of multi-mle to true
Multi_mle = Simulation_Results.MLE; % light mle vs time of multi-mle to true
%%
Uni_MLEs = cell(1,numel(Modules)); % will contain the uni-module mles
RM_times = []; % will contain the times to remove where we had 0 probability (no spikes emitted)
%% Getting speed
speed=data.task(L_D+1).tracking.MY_speed; % recorded speed in this experiment
speed_idc_below=find(speed<speed_thresh); % indices below speed threshold
%% Uni-module mles
cd(r_dir);
cd('Uni');

for mo=Modules
    name = [Na,sprintf(':module=%d,Tau=%dms,true_loc=0,max_shift=0,tr=1.mat',mo,tau)];
    load(name);       
    Uni_MLEs{mo} = Simulation_Results.MLE; % uni module mle
    
    if mo==1 % same Rel_Times for all modules
        Rel_times = Simulation_Results.Rel_times; % how much time backward did we look at?
        rm_times_multi_mle = Simulation_Results.rm_times_multi_mle; % times with no multi module mle
    end

    RM_times=[RM_times,Simulation_Results.rm_times]; % times with no uni mle for the mo'th module
end

RM_times=[RM_times,rm_times_multi_mle]; % times to be removed since no spikes were emitted in uni module decoding or no multi module mle
RM_times(RM_times==0)=[];

cd(ori_dir);
%% Removing irrelevant times
% Remove speed threshold, 0 prorbability from all modules (and when no multi module mle) and times before rel_times
Remove = unique([speed_idc_below',RM_times,1:Rel_times]); % all time points to be removed

for mo=Modules    
    M = Uni_MLEs{mo}; % uni-mle of the mo'th modules
    M(:,Remove) = []; % filtering
    Uni_MLEs{mo}=M; % putting back
end

Multi_err_to_true(Remove)=[]; % removing the same times from multi module error
Multi_mle(:,Remove)=[]; % removing the same times from multi module mle

%% Distances between uni-mle to multi-mle (and between uni pairs)
D_Uni_to_Mulit = zeros(numel(Modules),size(Multi_mle,2)); % will contain ditances between uni to multi after filtered times
D_Uni_to_Uni = zeros(pairs,size(Multi_mle,2)); % will contain ditances between uni pairs after filtered times

mo_Pairs = zeros(pairs,2); % will contain the pairs numbers
counter=1;

for mo=Modules
    
    u_mle = Uni_MLEs{mo}; % uni-mle of the mo'th modules
    D_Uni_to_Mulit(mo,:) = sqrt( (Multi_mle(1,:) - u_mle(1,:)).^2 + (Multi_mle(2,:) - u_mle(2,:)).^2 ); % distances for uni to multi


    for mo2=mo+1:numel(Modules) % running on second pair
        u_mle2 = Uni_MLEs{mo2}; % mle of second module in the pair
        
        D_Uni_to_Uni(counter,:) = sqrt ( (u_mle(1,:) - u_mle2(1,:)).^2  +  (u_mle(2,:) - u_mle2(2,:)).^2); % distance between mle's of a two modules vs time

        mo_Pairs(counter,:) = [mo,mo2]; % register pairs identities
        counter=counter+1;
    end

end
%% Histograms
figure(1);
H1 = histogram2(Multi_err_to_true,D_Uni_to_Mulit(1,:),'normalization','probability');
H1.BinWidth(1)=delta_err_multi;
V1 = H1.Values;
xlabel('Multi-module to true');
ylabel('Uni-module to Multi-mle');


figure(2);
H2 = histogram2(Multi_err_to_true,D_Uni_to_Mulit(2,:),'normalization','probability');
H2.BinWidth(1)=delta_err_multi;
V2 = H2.Values;
xlabel('Multi-module to true');
ylabel('Uni-module to Multi-mle');

if Rat~=3
    figure(3);
    H3 = histogram2(Multi_err_to_true,D_Uni_to_Mulit(3,:),'normalization','probability');
    H3.BinWidth(1)=delta_err_multi;    
    V3 = H3.Values;
    xlabel('Multi-module to true');
    ylabel('Uni-module to Multi-mle');
    
    S3 = sum(V3,2); % probability summed on MAE, 3
    I3 = find(cumsum(S3)>=thresh); % identical to I1
end

%% Plot conditional probability density for each module
S1 = sum(V1,2); % probability summed on MAE, 1
S2 = sum(V2,2); % probability summed on MAE, 2

I1 = find(cumsum(S1)>=thresh); % it's the same for all modules
I2 = find(cumsum(S2)>=thresh); % identical to I1

VV = {V1,V2};
SS = [S1,S2];

Uni_Bin = [H1.BinWidth(2), H2.BinWidth(2)]; % bin size along uni module axis
if Rat~=3 % if it's not Coconut so we hae 3 modules
    VV = {V1,V2,V3};
    SS = [S1,S2,S3];
    Uni_Bin = [H1.BinWidth(2), H2.BinWidth(2), H3.BinWidth(2)];
end

I = I1(1); % point where we stop
% delta_err = H1.BinWidth(1); % delta of error between 2 consecutive bins
% close all;

figure;
set(gcf, 'Position',[200, 300 1200 400])
Leg_Cell = cell(1,I); % will contain the legend
for m=1:numel(Modules)
    
    vv = VV{m};
    ss = SS(:,m);
    
    subplot(1,numel(Modules),m); hold all;
    for k=1:I % runnning on relevant multi-module err to plot
        plot(0:Uni_Bin(m):Uni_Bin(m)*(size(vv,2)-1) , vv(k,:)./max(vv(k,:))); % normalized conditional prbability
%         Leg_Cell{k} = ['MAE = ',num2str((k-1)*delta_err_multi),'-',num2str(k*delta_err_multi),'cm'];
        Leg_Cell{k} = [num2str((k-1)*delta_err_multi),'-',num2str(k*delta_err_multi),'cm'];
        xlim([0,Uni_Bin(m)*(size(vv,2)-1)]);
    end
    box on; grid on;
    xlabel(['\delta_',sprintf('%d [cm]',m)]);
%     xlabel([sprintf('u_%d',m), ' to multi distance [cm]']);
    ylabel('Normalized frequency');
    legend(Leg_Cell,'Location','bestoutside');
%     legend(Leg_Cell);
    set(gca,'fontsize',16);
       
end

%% Histograms for uni pairs
figure;
H1 = histogram2(Multi_err_to_true,D_Uni_to_Uni(1,:),'normalization','probability');
H1.BinWidth(1)=delta_err_multi;
V1 = H1.Values;
xlabel('Multi-module to true');
ylabel('Uni-module to Uni-module');

if Rat~=3 % Coconut has only a single pair
    figure;
    H2 = histogram2(Multi_err_to_true,D_Uni_to_Uni(2,:),'normalization','probability');
    H2.BinWidth(1)=delta_err_multi;    
    V2 = H2.Values;
    xlabel('Multi-module to true');
    ylabel('Uni-module to Uni-module');

    S2 = sum(V2,2); % probability summed on MAE, 2
    I2 = find(cumsum(S2)>=thresh); % identical to I1

    figure;
    H3 = histogram2(Multi_err_to_true,D_Uni_to_Uni(3,:),'normalization','probability');
    H3.BinWidth(1)=delta_err_multi;
    V3 = H3.Values;
    xlabel('Multi-module to true');
    ylabel('Uni-module to Uni-module');
    
    S3 = sum(V3,2); % probability summed on MAE, 3
    I3 = find(cumsum(S3)>=thresh); % identical to I1
end

%% Plot conditional probability density for each module
S1 = sum(V1,2); % probability summed on MAE, 1
I1 = find(cumsum(S1)>=thresh); % it's the same for all modules

VV = {V1};
SS = S1;

Uni_Bin = [H1.BinWidth(2), H2.BinWidth(2)]; % bin size along uni module axis
if Rat~=3 % if it's not Coconut so we has 3 modules
    VV = {V1,V2,V3};
    SS = [S1,S2,S3];
    Uni_Bin = [H1.BinWidth(2), H2.BinWidth(2), H3.BinWidth(2)];
end

I = I1(1); % point where we stop
% delta_err = H1.BinWidth(1); % delta of error between 2 consecutive bins
% close all;

figure;
set(gcf, 'Position',[200, 300 1200 400])
Leg_Cell = cell(1,I); % will contain the legend
for m=1:pairs % running on all pairs
    
    vv = VV{m};
    ss = SS(:,m);
    
    subplot(1,pairs,m); hold all;
    for k=1:I % runnning on relevant multi-module err to plot
        plot(0:Uni_Bin(m):Uni_Bin(m)*(size(vv,2)-1),vv(k,:)./max(vv(k,:))); % normalized conditional prbability
%         Leg_Cell{k} = ['MAE = ',num2str((k-1)*delta_err_multi),'-',num2str(k*delta_err_multi),'cm'];
        Leg_Cell{k} = [num2str((k-1)*delta_err_multi),'-',num2str(k*delta_err_multi),'cm'];
        xlim([0,Uni_Bin(m)*(size(vv,2)-1)]);
    end
    box on; grid on;
%     xlabel(['u_',num2str(mo_Pairs(m,1)),' to u_',num2str(mo_Pairs(m,2)),' distance [cm]']);
    xlabel(['\Delta_',sprintf('%d_{,%d} [cm]',mo_pairs(m,1),mo_pairs(m,2))]);
    ylabel('Normalized frequency');
    legend(Leg_Cell);
%     legend(Leg_Cell,'Location','bestoutside');
    set(gca,'fontsize',16);
       
end