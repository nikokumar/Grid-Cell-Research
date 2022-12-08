clear; close all;
Rat=2; % 1 for Bubble, 11for Bubble 2nd dataset, 2 for Roger,3 for Coconut, 4 for Rolf
L_D=0; % 0 for dark and 1 for light
speed_thresh=3; % speed threshold [cm/sec]

tau=100; % the used exp time constant
Modules=1:3;

pairs=factorial(numel(Modules))/2; % # of pairs of modules
tr=1; % the trial we look at, they are identical so it doesn't matter

%% Errors and segment length
seg_min_time = 1; % segment minimum continuous time [sec]

% Higher bound of above error will cover >=80-85% of histogram probability density so have differnt maximal value for light and dark
if Rat==1 % Bubble
    e_thresh_below = 10; % maximal value for low err segments [cm]
    if L_D==0 % dark
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 51; % maximal value for err for high err segments [cm] 57cm>=85%, 51cm>=80%
    elseif L_D==1 % light
        e_thresh_above_min = 15; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 30; % maximal value for err for high err segments [cm] 30cm=85%
    end
elseif Rat==11 % Bubble 2nd dataset
    e_thresh_below = 10; % maximal value for low err segments [cm]
    if L_D==0 % dark
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 48; % maximal value for err for high err segments [cm] 54cm>=85%, 48cm>=80%
    elseif L_D==1 % light
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 21; % maximal value for err for high err segments [cm] 21cm>=85%
    end
elseif Rat==2 % Roger
    e_thresh_below = 10; % maximal value for low err segments [cm]
    if L_D==0 % dark
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 45; % maximal value for err for high err segments [cm] 72cm>=85%, 45cm>=80%
    elseif L_D==1 % light
        e_thresh_above_min = 15; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 21; % maximal value for err for high err segments [cm] 21cm>=85%
    end
elseif Rat==3 % Coconut
    e_thresh_below = 10; % maximal value for low err segments [cm]
    if L_D==0 % dark
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 102; % maximal value for err for high err segments [cm] 105cm>=85%, 102cm>=80%
    elseif L_D==1 % light
        e_thresh_above_min = 15; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 60; % maximal value for err for high err segments [cm]
    end
elseif Rat==4
    e_thresh_below = 10; % maximal value for low err segments [cm]
    if L_D==0 % dark
        e_thresh_above_min = 20; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 66; % maximal value for err for high err segments [cm] 66cm>=80%, 75cm>=85%
    elseif L_D==1 % light
        e_thresh_above_min = 15; % minimal value for err for high err segments [cm]
        e_thresh_above_max = 33; % maximal value for err for high err segments [cm]
    end
end
% e_thresh_above_max = 45; % same maximal value for err for high err segments for light and dak[cm]

ori_dir=pwd;
%%
if L_D==0
    Na='Dark';
elseif L_D==1
    Na='Light';
end

if Rat==1 % Bubble
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data');
    ra='Bubble';
elseif Rat==11 % Bubble 2nd dataset
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Bubble_data 2');
    ra='Bubble 2';
elseif Rat==2 % Roger
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Roger_data');
    ra='Roger';
elseif Rat==3 % Coconut
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Coconut_data')
    ra='Coconut';
    Modules=1:2;
elseif Rat==4 % Rolf
    cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Rat_Rolf_data')
    ra='Rolf';
end
%% Loading speed
load('data.mat'); % loading the relevant data file
speed=data.task(L_D+1).tracking.MY_speed; % recorded speed
speed_idc_above=find(speed>=speed_thresh); % indices above speed threshold

Time = data.task(L_D+1).tracking.timestamp; % corresponding time [sec]
dt = mean(Time(2:end)-Time(1:end-1)); % time interval between 2 time points [sec]

seg_min_time = round(seg_min_time/dt)-1;
if seg_min_time<=0 % if we use minimal segment time
    seg_min_time=1;
end
cd(ori_dir);
%% Smoothing parameter
sigma = .05; % desired sigma of gaussian in sec [50msec]
sigma = round(sigma/dt); % sigma in indices (dt units)

fac=5; % factor between std and window size of smoothdata gaussian function
win = fac*sigma; % the window given as an input to the smoothdata function

%% Loading  uni-module error
cd('/Users/haggaiagmon/Documents/MATLAB/Torgeir data New/Exp decoder');
cd(ra);
cd('0 shifts results for focus');
r_dir=pwd;
cd('Uni');

Uni_MLEs=cell(1,numel(Modules)); % will contain the mle for each module
RM_times=[]; % will contain the times to remove where we had 0 probability (no spikes emitted)

for mo=Modules
    name = [Na,sprintf(':module=%d,Tau=%dms,true_loc=0,max_shift=0,tr=%d.mat',mo,tau,tr)]; % all trials are identical
    load(name);

    Uni_MLEs{1,mo} = Simulation_Results.MLE; % estimated uni position

    RM_times=[RM_times,Simulation_Results.rm_times]; % times with no uni mle for the mo'th module
    if mo==1 % if unit cell was around multi module mle there could be also times where there was no multi-mle, identical for all modules
        rm_times_multi_mle = Simulation_Results.rm_times_multi_mle; % times with no multi module mle
        Rel_times= Simulation_Results.Rel_times; % times to crop from start
    end
end

RM_times=[RM_times,rm_times_multi_mle]; % times to be removed since no spikes were emitted in uni module decoding or no multi module mle
RM_times(RM_times==0)=[];

%% Multi-module MLE
cd(r_dir);
cd('Multi');

name = [Na,sprintf(':Tau=%dms,spatial_shift=0,tr=%d.mat',tau,tr)]; % all trials are identical
load(name);
multi_mle = Simulation_Results.MLE;

%% Error and smoothing of multi-module 
multi_err = Simulation_Results.ERR; % mean absoulute error
[multi_err,~] = smoothdata(multi_err,2,'gaussian',win); % smoothing
%% Times with error above/below the error threshold 
error_idc_above_min = find(multi_err >= e_thresh_above_min); % times with error above minimal value for err for high err segments
error_idc_above_max = find(multi_err < e_thresh_above_max); % times with error below maximal value for err for high err segments
error_idc_above = intersect(error_idc_above_min,error_idc_above_max); % times with error higher e_thresh_above_min and lower than e_thresh_above_max


error_idc_below = find(multi_err <= e_thresh_below); % times with error below min error
%% Relevant times
% Times which decoding didn't skipped, and had enough speed
valid = 1:numel(multi_err);
valid([1:Rel_times,RM_times])=[]; % this are valid times with decoding result

Inter = intersect(speed_idc_above',valid); % time points which are above minimum speed and wern't skipped

Times_above_err = intersect(Inter,error_idc_above); % intersecting with time points above minimal error
Times_below_err = intersect(Inter,error_idc_below); % intersecting with time points below minimal error
%% Segments
contin_above = Times_above_err(seg_min_time+1:end) - Times_above_err(1:end-seg_min_time); % is it conitouous ?
seg_start_above = find(contin_above==seg_min_time); % starting time of all segments above error threshold

Seg_Start_Above = zeros(1,numel(seg_start_above)); % will contain stranger segments starting times
counter=1;
flag=1; % flags if segment is to be used or not
for i=1:numel(seg_start_above)-1
    if flag==1 % if its a segment to save
        if i==1
            Seg_Start_Above(counter) = seg_start_above(i); % saving initial segment
        elseif i~=1
            Seg_Start_Above(counter) = seg_start_above(i-1); % saving segments
        end
        flag=0; % our defult now is not to save the next segment
        counter=counter+1;
    elseif flag==0 % we need to check if it's a segment to be saved
        d_seg = seg_start_above(i) - Seg_Start_Above(counter-1); % difference between segments starting time
        if d_seg>=seg_min_time % only if if it's more the minimal time (non-overlapping time) we flag it
            flag=1;
        end
    end
end
Seg_Start_Above(Seg_Start_Above==0)=[];
n_seg_above=numel(Seg_Start_Above); % # of segments above error threshold that stand in this conditions



contin_below = Times_below_err(seg_min_time+1:end) - Times_below_err(1:end-seg_min_time); % is it conitouous ?
seg_start_below = find(contin_below==seg_min_time); % starting time of all segments below error threshold

Seg_Start_Below = zeros(1,numel(seg_start_below)); % will contain stranger segments starting times
counter=1;
flag=1; % flags if segment is to be used or not
for i=1:numel(seg_start_below)-1
    if flag==1 % if its a segment to save
        if i==1
            Seg_Start_Below(counter) = seg_start_below(i); % saving initial segment
        elseif i~=1
            Seg_Start_Below(counter) = seg_start_below(i-1); % saving segments
        end
        flag=0; % our defult now is not to save the next segment
        counter=counter+1;
    elseif flag==0 % we need to check if it's a segment to be saved
        d_seg = seg_start_below(i) - Seg_Start_Below(counter-1); % difference between segments starting time
        if d_seg>=seg_min_time % only if if it's more the minimal time (non-overlapping time) we flag it
            flag=1;
        end
    end
end
Seg_Start_Below(Seg_Start_Below==0)=[];    
n_seg_below=numel(Seg_Start_Below); % # of segments below error threshold that stand in this conditions

%% Distances between unis to multi and between unis vs time    
d_uni_to_multi = zeros(numel(Modules),size(multi_mle,2)); % will contain distance between uni to multi vs time
d_between_unis = zeros(pairs,numel(multi_err)); % will contain distance between unis vs time

counter=1;
for m1=1:numel(Modules) % running on first pair

    mle1 = Uni_MLEs{m1}; % mle of first module in the pair

    d_uni_to_multi(m1,:) = sqrt ( (mle1(1,:) - multi_mle(1,:)).^2  +  (mle1(2,:) - multi_mle(2,:)).^2); % uni-to multi difference vs time

    for m2=m1+1:numel(Modules) % running on second pair
        mle2 = Uni_MLEs{m2}; % mle of second module in the pair

        d_between_unis(counter,:) = sqrt ( (mle1(1,:) - mle2(1,:)).^2  +  (mle1(2,:) - mle2(2,:)).^2); % mean distance between mle's or a two modules for this segment

        counter=counter+1;
    end
end

%% Mean distnaces in the given segments
Distances_uni_to_multi_above = zeros(numel(Modules),n_seg_above); % will contain the mean distances between each uni module mle to multi module mle for segments above error threshold
distances_between_unis_above = zeros(pairs,n_seg_above); % will contain the mean distnace between any pair of modules mles' vs time for segments above error threshold
distance_multi_to_true_above = zeros(1,n_seg_above); % will contain the mean distances between multi module mle to true location for segments above error threshold



Distances_uni_to_multi_below = zeros(numel(Modules),n_seg_below); % will contain the mean distances between each uni module mle to multi module mle for segments below error threshold
distances_between_unis_below = zeros(pairs,n_seg_below); % will contain the mean distnace between any pair of modules mles' vs time for segments below error threshold
distance_multi_to_true_below = zeros(1,n_seg_below); % will contain the mean distances between multi module mle to true location for segments below error threshold


for k=1:n_seg_above % running on all segments above err thresohld
    start_point_above = Seg_Start_Above(k); % start point of segment    

    for j=1:size(Distances_uni_to_multi_above,1) % run on uni to multi
        Distances_uni_to_multi_above(j,k) = mean(d_uni_to_multi(j,Times_above_err(start_point_above):Times_above_err(start_point_above)+seg_min_time));
    end
    for j=1:size(distances_between_unis_above,1) % run on uni to uni
        distances_between_unis_above(j,k) = mean(d_between_unis(j,Times_above_err(start_point_above):Times_above_err(start_point_above)+seg_min_time));
    end
    
    distance_multi_to_true_above(k) = mean(multi_err(Times_above_err(start_point_above):Times_above_err(start_point_above)+seg_min_time));
end

for k=1:n_seg_below % running on all segments below err thresohld
    start_point_below = Seg_Start_Below(k); % start point of segment    

    for j=1:size(Distances_uni_to_multi_below,1) % run on uni to multi
        Distances_uni_to_multi_below(j,k) = mean(d_uni_to_multi(j,Times_below_err(start_point_below):Times_below_err(start_point_below)+seg_min_time));
    end
    for j=1:size(distances_between_unis_below,1) % run on uni to uni
        distances_between_unis_below(j,k) = mean(d_between_unis(j,Times_below_err(start_point_below):Times_below_err(start_point_below)+seg_min_time));
    end
    
    distance_multi_to_true_below(k) = mean(multi_err(Times_below_err(start_point_below):Times_below_err(start_point_below)+seg_min_time));
end
%% Results from original dataset
cd(ori_dir);
disp(['Distances unis to multi Above: ',Na]);
disp(mean(Distances_uni_to_multi_above,2));
disp('# of segments Above:');
disp(n_seg_above);
disp('Multi-mle to true distance Above:');
disp(mean(distance_multi_to_true_above));

disp(['Distances unis to multi Below: ',Na]);
disp(mean(Distances_uni_to_multi_below,2));
disp('# of segments Below:');
disp(n_seg_below);
disp('Multi-mle to true distance Below:');
disp(mean(distance_multi_to_true_below));


disp([Na,' overall mean error between multi-mle to true location:']);
disp( (n_seg_above*mean(distance_multi_to_true_above) + n_seg_below*mean(distance_multi_to_true_below))/(n_seg_above+n_seg_below));

%% Results Uni to Multi
Vals_U_to_M_above = mean(Distances_uni_to_multi_above,2); % mean values for unis to multi above
Err_U_to_M_above = std(Distances_uni_to_multi_above,[],2)/sqrt(size(Distances_uni_to_multi_above,2)); % corresponding error (sem)
Vals_U_to_M_below = mean(Distances_uni_to_multi_below,2); % mean values  for unis to multi below
Err_U_to_M_below = std(Distances_uni_to_multi_below,[],2)/sqrt(size(Distances_uni_to_multi_below,2)); % corresponding error (sem)

Vals_M_to_T_above = mean(distance_multi_to_true_above); % mean values for multi to true above
Err_M_to_T_above = std(distance_multi_to_true_above)/sqrt(numel(distance_multi_to_true_above));  % corresponding error (sem)
Vals_M_to_T_below = mean(distance_multi_to_true_below); % mean values for multi to true below
Err_M_to_T_below = std(distance_multi_to_true_below)/sqrt(numel(distance_multi_to_true_below));  % corresponding error (sem)

Overall_err = (n_seg_above*Vals_M_to_T_above + n_seg_below*Vals_M_to_T_below)/(n_seg_above+n_seg_below);
%% Plot bars Uni to Multi
br_widt=.8;

X_name = categorical({'\delta_1','\delta_2','\delta_3','MAE'});
X_name = reordercats(X_name,{'\delta_1','\delta_2','\delta_3','MAE'});

if Rat==3
    X_name = categorical({'\delta_1','\delta_2','MAE'});
    X_name = reordercats(X_name,{'\delta_1','\delta_2','MAE'});
end


figure(1); 
b = bar(X_name,[Vals_U_to_M_below',Vals_M_to_T_below ; Vals_U_to_M_above',Vals_M_to_T_above],br_widt);
X_pos=zeros(2,4);
if Rat==3
    X_pos=zeros(2,3);
end
for a_b=1:2 % get x locations of bars for aboe abd below
    X_pos(a_b,:) = b(a_b).XEndPoints; % x position of bars    
end
% X_pos(1,:)=X_pos(1,:)-1;
% X_pos(2,:)=X_pos(2,:)+1;

hold on;
errorbar(X_pos,[Vals_U_to_M_below',Vals_M_to_T_below ; Vals_U_to_M_above',Vals_M_to_T_above],[Err_U_to_M_below',Err_M_to_T_below ; Err_U_to_M_above',Err_M_to_T_above],'k.')
ylabel('Distance [cm]');
legend('Small error periods','Large error periods','Location','northwest');
title(Na);
set(gca,'fontsize',16); grid on;

%% Results Uni to Uni
Vals_U_to_U_above = mean(distances_between_unis_above,2); % mean values for unis to multi above
Err_U_to_U_above = std(distances_between_unis_above,[],2)/sqrt(size(distances_between_unis_above,2)); % corresponding error (sem)
Vals_U_to_U_below = mean(distances_between_unis_below,2); % mean values  for unis to multi below
Err_U_to_U_below = std(distances_between_unis_below,[],2)/sqrt(size(distances_between_unis_below,2)); % corresponding error (sem)

%% Plot bars Uni to Uni
br_widt=.8;

X_name = categorical({'\Delta_{1,2}','\Delta_{1,3}','\Delta_{2,3}','MAE'});
X_name = reordercats(X_name,{'\Delta_{1,2}','\Delta_{1,3}','\Delta_{2,3}','MAE'});

if Rat==3
    X_name = categorical({'\Delta_{1,2}','MAE'});
    X_name = reordercats(X_name,{'\Delta_{1,2}','MAE'});
end


figure(2); 
b = bar(X_name,[Vals_U_to_U_below',Vals_M_to_T_below ; Vals_U_to_U_above',Vals_M_to_T_above],br_widt);
X_pos=zeros(2,4);
if Rat==3
    X_pos=zeros(2,3);
end
for a_b=1:2 % get x locations of bars for aboe abd below
    X_pos(a_b,:) = b(a_b).XEndPoints; % x position of bars    
end
% X_pos(1,:)=X_pos(1,:)-1;
% X_pos(2,:)=X_pos(2,:)+1;

hold on;
errorbar(X_pos,[Vals_U_to_U_below',Vals_M_to_T_below ; Vals_U_to_U_above',Vals_M_to_T_above],[Err_U_to_U_below',Err_M_to_T_below ; Err_U_to_U_above',Err_M_to_T_above],'k.')
ylabel('Distance [cm]');
legend('Small error periods','Large error periods','Location','northwest');
title(Na);
set(gca,'fontsize',16); grid on;