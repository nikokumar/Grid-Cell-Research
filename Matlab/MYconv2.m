function [ conv2d ] = MYconv2( A,B )
% A and B are 2d squared matrices which will get convolved using 2d fft
% Yields the 'same' size convolution for 2 square matrices
rowA=size(A,1); colA=size(A,2); % rowA = colA
rowB=size(B,1); colB=size(B,2); % rowB = colB;

A1 = [A zeros(rowA,colB-1); zeros(rowB-1,colA+colB-1)];

B1 = zeros(size(A1));    B1(1:rowB,1:colB) = B;

conv2dfull = ifft2(fft2(A1).*fft2(B1)); % the convolution with fft, full


conv2dfull(conv2dfull<0)=0; % ifft(fft) produces numerical errors, negative numbers


c1=size(conv2dfull,1);  % = a1+b1-1, always odd

% Con = conv2(A,B,'same'); % the convolution with matlab conv2 function, same

%% making the fft conv from full to same
if mod(rowA,2)==0
    cut = ceil((c1-rowA)/2);
    conv2d = conv2dfull(cut+1:cut+rowA,cut+1:cut+colA);
else
    cut=(c1-rowA)/2;
    conv2d = conv2dfull(cut+1:cut+rowA,cut+1:cut+colA);
end
%% making sure we get the correct result
% disp('convolution with matlab CONV2:');
% disp(Con);
% disp('conv2dfull');
% disp(conv2dfull);
% disp('convolution with my FFT:');
% disp(conv2d);
end