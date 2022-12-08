function [RF,dr] = GaussianRF(lambda,L,theta )
%% Two approaches to get different phases
% first appraoch: create grid twice larger and then cut in halves to get
% different phases, can rotate the grid before

% rangemax=2*rangemax; % we do areana 2 times bigger because in MM_GaussianRF we will cut it in half to generate different phases 
% rangemin=-rangemax; %[meters]
% deltarange=0.01; % resoultion [meters]
% [X,Y]=meshgrid(rangemin:deltarange:rangemax,rangemin:deltarange:rangemax); % the arena
% theta= rand.*2.*pi; % random angle to rotate the grid
% RF=zeros(numel(rangemin:deltarange:rangemax),numel(rangemin:deltarange:rangemax));


% Second approach: do not rotate grid (theta=0) and elongate x axis by
% lambda and y axis by sqrt(3)*lambda and crop randomly accordingly

minus_L=-L; % [cm]
dr=1; %0.0025 resolution [cm]

[X,Y]=meshgrid(minus_L:dr:L+lambda,minus_L:dr:L+(sqrt(3)/2)*lambda); % the arena
% theta=pi/2;%0;
RF=zeros(numel(minus_L:dr:L+(sqrt(3)/2)*lambda),numel(minus_L:dr:L+lambda));

%% Creating the centers of the gaussians
% Xcenter1 will go with Ycenter1, coordinates of centers of gaussians
% Xcenter2 will go with Ycenter2, coordinates of centers of gaussias

% making the Xcenters
v=0;
i=1;
xplus=zeros(ceil(L/lambda)+2,1); % including the origin
while v<L+2*lambda % add anther lambda because we will rotate the RF
    xplus(i,1)=v;
    v=v+lambda;
    i=i+1;
end

v=-lambda;
i=1;
xminus=zeros(ceil(abs(minus_L)/lambda)+1,1); % not including the origin
while v>minus_L-2*lambda % add anther lambda because we will rotate the RF
    xminus(i,1)=v;
    v=v-lambda;
    i=i+1;
end
xminus(xminus==0)=[]; % excluding the origin

Xcenter1=[flipud(xminus);xplus];

% need to make Xcenter2 formally and not from Xcenter1
v=lambda/2;
i=1;
xplus2=zeros(ceil(L/lambda)+1,1);
while v<L+2*lambda % add anther lambda because we will rotate the RF
    xplus2(i,1)=v;
    v=v+lambda;
    i=i+1;
end

v=-lambda/2;
i=1;
xminus2=zeros(ceil(abs(minus_L)/lambda)+1,1);
while v>minus_L-2*lambda % add anther lambda because we will rotate the RF
    xminus2(i,1)=v;
    v=v-lambda;
    i=i+1;
end

Xcenter2=[flipud(xminus2);xplus2];
Xcenter2(Xcenter2==0)=[]; % excluding mistake of the origin

% making the Ycenters
a=2*sin(pi./3)*lambda;
v=0;
i=1;
yplus=zeros(ceil(L/a)+1,1); % including the origin
while v<L+2*lambda % add anther lambda because we will rotate the RF
    yplus(i,1) = v;
    v=v+a;
    i=i+1;
end
v=-a;
i=1;
yminus=zeros(ceil(abs(minus_L)/a)+1,1); % not including the origin
while v>minus_L-2*lambda % add anther lambda because we will rotate the RF
    yminus(i,1)=v;
    v=v-a;
    i=i+1;
end
yminus(yminus==0)=[]; % excluding origin

Ycenter1=[flipud(yminus);yplus];

% need to make Ycenter2 formally and not from Ycenter1
v=a/2;
i=1;
yplus2=zeros(ceil(L/a)+1,1);
while v<L+2*lambda % add anther lambda because we will rotate the RF
    yplus2(i,1)=v;
    v=v+a;
    i=i+1;
end

v=-a/2;
i=1;
yminus2=zeros(ceil(abs(minus_L)/a)+1,1);
while v>minus_L-2*lambda % add anther lambda because we will rotate the RF
    yminus2(i,1)=v;
    v=v-a;
    i=i+1;
end
Ycenter2= [flipud(yminus2);yplus2];
Ycenter2(Ycenter2==0)=[]; % excluding mistake of the origin
%% Mu1 and Mu2 are vectors containing the centers of all gaussians in the grid
ye1=numel(Ycenter1);
xe1=numel(Xcenter1);
Mu1=zeros(xe1.*ye1,2);

v=1;
u=ye1;
for i=1:xe1
    Mu1(v:u,1)=Xcenter1(i);
    Mu1(v:u,2)=Ycenter1;
    v=v+ye1;
    u=u+ye1;
end

ye2=numel(Ycenter2);
xe2=numel(Xcenter2);
Mu2=zeros(xe2.*ye2,2);

v=1;
u=ye2;
for i=1:xe2
    Mu2(v:u,1)=Xcenter2(i);
    Mu2(v:u,2)=Ycenter2;
    v=v+ye2;
    u=u+ye2;
end
%% Generating gaussians with the given centers

sigmaX = 0.015*(lambda^2); sigmaY =  0.015*(lambda^2); % the variances of the gaussians are related to lambda, [lambda]=meters
Sigma=[sigmaX, 0; 0, sigmaY]; % Covariance matrix
MU=[Mu1;Mu2]; % vector containing all the centers of all of the gaussians in the grid

for i=1:size(MU,1)
    
    gaussian = mvnpdf([cos(theta).*X(:)-sin(theta).*Y(:) sin(theta).*X(:)+cos(theta).*Y(:)],MU(i,:),Sigma);
    gaussian = reshape(gaussian, size(X,1), size(X,2));
   
    RF=RF+gaussian;
end
% figure
% imagesc(RF);
end