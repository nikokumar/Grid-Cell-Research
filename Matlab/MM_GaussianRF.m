function [ RFS,dr,RND,L,Lambdas] = MM_GaussianRF( Modules, N, L, Rat )
% this code will generate the receptive fields of all modules according to gaussian receptive field.
% then it will create different phases for all different neurons in this moudle

% Inputs: Moudle = number of different modules
% n = number of neurons in each module
% L = radi of the arena [cm]

% Output: RFS = cell, each element in this cell contain the recptive field of a specific module
% Rnd = vector containing the uniform phases, coordinates of uniform phases
% on the receptive field in RFS
% L=[Ly,Lx] length of small arean in Y an X


RFS=cell(1,Modules); % will contain the general receptive field for each module

%% Generating the grid spacing in the vector Lambdas
% a=sqrt(2); % ratio between module's spacings'
% spacing_const=25; % 
% Lambdas = a.^(1:Modules).*spacing_const; % vector containing the different grid spacing [cm]
% the first element is the smallest lambda

if Rat==1 || Rat==11
    Lambdas =[45,65,95]; % mimic Bubble spacings
elseif Rat==2
    Lambdas =[50,80,100]; % mimic Roger spacings
end
%% Generating the gaussian receptive fields for each moudle + uniform phases
RND=cell(1,Modules); % will contain phases for all modules
for i=1:Modules % running on all modules
    
    lambda=Lambdas(i); % look at the given i'th lambda
    
%     rng(i);
    theta=rand*pi; % random orientation angle for the i'th module
    
    [RF,dr] = GaussianRF( lambda,L,theta ); % RF is a "big" receptive field for the given lambda
    
    RF=RF./(max(max(RF))); % making the maximum of the recptive field = 1
    
    RFS{i}=RF;   

    %% Generating unifrom phases (randomly ?)
    % putting directly maximum # of neurons that we can uniformly, and the rest randomly
    PixelsX=floor(lambda/dr);               % # of pixels in X
    PixelsY=floor((lambda*(sqrt(3)/2))/dr); % # of pixels in Y
    %% Complete random phases
    Rnd=zeros(N(i),2); % Rnd(i,:) coordinates of start point of small grid for the i'th neuron
    Rnd(:,1)=ceil(PixelsX.*rand(N(i),1));
    Rnd(:,2)=ceil(PixelsY.*rand(N(i),1));
    RND{i}=Rnd;
end

end