function [TCs,STs] = Sythetic_TC_and_poiss_spikes(Rat,X_traj,Y_traj,Modules,N,dt,Lim, identical,max_shift,tr)
% Generate synthetic spike train as a function of arena size.
% It is used for measureing the cross correlation between inter-module neurons:
% correlation should decrease as arena grows

rng(tr); % random generator for poisson spikes
Time = numel(X_traj);
max_fr=30; % maxial firing rate of all neurons [Hz]
%% Generate tuning curves
% Receptive fileds, spatial resulotion and phases for the given # of modules, allocation of neurons and arena size
[ RFS,dr,RND,~,~] = MM_GaussianRF( Modules, N, Lim, Rat );

[X,Y] = meshgrid(-Lim:dr:Lim, -Lim:dr:Lim);
R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>Lim^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

%% Generate spike trains for this trajectory
STs=cell(1,Modules); % will contain the spike trains for each module
TCs=cell(1,Modules); % will contain the corrsponding tuning curves

Xwalkround=round(X_traj); % rounding to cm resolution, x coordinate
Ywalkround=round(Y_traj); % rounding to cm resolution, y coordinate
MapXwalk = Xwalkround + Lim; % X index to use on the RF matrix
MapYwalk = Ywalkround + Lim; % Y index to use on the RF matrix

%% Generate spatial shifts for spike generation
shifts_x = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in x axis
shifts_y = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in y axis

if identical==1 % if we actually perform identical spatial shifts
    rng(tr*1000); % new seed for identical shifts
    shifts_x = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in x axis
    shifts_y = round(2.*max_shift.*rand(1,3) - max_shift); % spatial shifts in y axis

    shifts_x(2:3) = shifts_x(1); % same x shifts for all
    shifts_y(2:3) = shifts_y(1); % same y shifts for all
end

% Updating phases
RND_Shifted = cell(1,Modules); % will contain shifted phases for all neurons in each module
for i = 1:Modules % running on all modules
    Rnd_shifted = RND{i}; % relevant phases
    Rnd_shifted(:,1) = Rnd_shifted(:,1) + shifts_x(i);
    Rnd_shifted(:,2) = Rnd_shifted(:,2) + shifts_y(i);

    RND_Shifted{i} = Rnd_shifted;
end
%%
for i=1:Modules % running on all modules
    
    Rnd_shifted=RND_Shifted{i}; % relevant shifted phases for generating spikes
    Rnd_ori=RND{i}; % relevant original phases for decoding

    RF=RFS{i}; % the BIG receptive field of this module, its normalized and so the max value=1
    st=zeros(N(i),Time); % will contain the spike trains vs time for the i'th module
    tc=zeros(2*Lim+1,2*Lim+1,N(i));
    
    for j=1:N(i) % running on all neurons in that module
        rf_ori = RF(Rnd_ori(j,2):Rnd_ori(j,2)+2*Lim,Rnd_ori(j,1):Rnd_ori(j,1)+2*Lim); % original receptive field of the j'th neuron for decoding
        rf_shifted = RF(Rnd_shifted(j,2):Rnd_shifted(j,2)+2*Lim,Rnd_shifted(j,1):Rnd_shifted(j,1)+2*Lim); % shifted receptive field of the j'th neuron for generating spikes
        % coordinates of Rnd=Rnd(x,y) so rows=y and colums=x
        rf_ori(idc)=0; % 0 firing rate outside the arena (the trajectory can't also exceed the wrena boundaries)
        rf_shifted(idc)=0; % 0 firing rate outside the arena (the trajectory can't also exceed the wrena boundaries)
        
        rf_shifted = rf_shifted.*max_fr; % setting maximal firing rate for rate maps for encoding spikes
        rf_ori = rf_ori.*max_fr; % setting maximal firing rate for rate maps for decoding spikes

        tc(:,:,j) = rf_ori; % register original tuning curve
        
        idx = sub2ind(size(rf_shifted),MapYwalk,MapXwalk); % translate trajectory long indices to short index
        fr_vs_time = rf_shifted(idx); % firing rate vs time of shifted tuning curves
        
        st(j,:) = poissrnd(fr_vs_time.*dt);  % poisson spikes      
    end
    STs{i}=st; % spike trains of all neurons from the i'th module
    TCs{i}=tc; % corresponding tuning curves from all neurons in the i'th module
end

end