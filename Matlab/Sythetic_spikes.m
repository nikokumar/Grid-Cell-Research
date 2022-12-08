function [STs,TCs] = Sythetic_spikes(Rat,Modules,N,L_D,dt,L,tr,cluster)
% Generate synthetic spike train as a function of arena size.
% It is used for measureing the cross correlation between inter-module neurons:
% correlation should decrease as arena grows
%% Hyper-parameters
% Using same speeds as in actual recording
% Rat=1; % 1 for Bubble, 2 for Roger
% L_D=0; % 0 for dark, 1 for light

% L = 50; % arena radi [cm]

rng(tr); % random generator

Time = 10^6;
% Time = 50*L^2; % # of time points, need more time points as arena grows

% Modules=3; % total # of modules
% N=[66,61,38]; % # of neurons in each module
max_fr=30; % maxial firing rate of all neurons [Hz]
%% Generate tuning curves
% Receptive fileds, spatial resulotion and phases for the gieven # of modules, allocation of neurons and arena size
[ RFS,dr,RND,~,~] = MM_GaussianRF( Modules, N, L, Rat );

[X,Y] = meshgrid(-L:dr:L, -L:dr:L);
R_sq=X.^2 + Y.^2; % squared distance from origin
[out_y,out_x]=find(R_sq>L^2); % coordinates to set to 0 firing rate since they are outsile of circle
idc = sub2ind(size(X), out_y,out_x); % linear indices of points outside of circled arena

%% Generate trajectory
[ X_R_traj, Y_R_traj] = RandomWalk(Rat,L_D,L,Time,cluster);
%% Generate spike trains for this trajectory
STs=cell(1,Modules); % will contain the spike trains for each module
TCs=cell(1,Modules); % will contain the corrsponding tuning curves

Xwalkround=round(X_R_traj); % rounding to cm resolution, x coordinate
Ywalkround=round(Y_R_traj); % rounding to cm resolution, y coordinate
MapXwalk = ((Xwalkround+L)./dr) + 1; % X index to use on the RF matrix
MapYwalk = ((Ywalkround+L)./dr) + 1; % Y index to use on the RF matrix


for i=1:Modules % running on all modules
    
    Rnd=RND{i}; % relevant phases
    RF=RFS{i}; % the BIG receptive field of this module, its normalized and so the max value=1
    st=zeros(N(i),Time); % will cotain the spike trains vs time for the i'th module
    tc=zeros(2*L+1,2*L+1,N(i));
    
    for j=1:N(i) % running on all neurons in that module
        rf = RF(Rnd(j,2):Rnd(j,2)+2*L,Rnd(j,1):Rnd(j,1)+2*L); % receptive field of the j'th neuron
        % coordinates of Rnd=Rnd(x,y) so rows=y and colums=x
        rf(idc)=0; % 0 firing rate outside the arena (the trajectory can't also exceed the wrena boundaries)
        
        tc(:,:,j) = rf;
        
        idx = sub2ind(size(rf),MapYwalk,MapXwalk); % translate trajectory long indices to short index
        fr_vs_time = rf(idx).*max_fr; % firing rate vs time
        
        st(j,:) = poissrnd(fr_vs_time.*dt);  % poisson spikes      
    end
    STs{i}=st; % spike trains of all neurons from the i'th module
    TCs{i}=tc; % corresponding tuning curves from all neurons in the i'th module
end

end