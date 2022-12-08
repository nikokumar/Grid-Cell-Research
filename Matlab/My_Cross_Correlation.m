function [C,lags] = My_Cross_Correlation(a,b,l,mat_norm)
% inputs: a and b are two vectors from the same dimension to cross correlate.
%         l is an integer = maximum of coordinate lags to compute.
%         mat_norm: =1 for matlab normalization [0,1], ~=1 for pearson correlation [-1,1]
% outputs: C is cross correlation, and lags are the lags...

lags=-abs(round(l)):1:abs(round(l)); % the lags
a=a(:)'; % row vector
b=b(:)'; % row vector

if l>=numel(a) % we can shift the vectros by more of their dimension
    l=numel(a)-1;
end

C=zeros(1,numel(lags)); % will contain the cross-correlation
counter=1;
%% Matlab normalization for correlation [0,1]
if mat_norm==1 % matalb normalization, max value=1, min value=0
    for k=-l:1:l % more than the zero time lag is requested

        if k<0 % negative time lag
            C(counter) = sum(a(abs(k)+1:end) .* b(1:end-abs(k)));
            counter=counter+1;
        end
        if k==0 % 0 time lag
            C(counter) = sum(a.*b);
            counter=counter+1;
        end
        if k>0 % positive time lag
            C(counter) = sum(a(1:end-k).*b(k+1:end));
            counter=counter+1;
        end
    end

    matlab_normalization=sqrt( sum(a.*a) .* sum(b.*b));
    C=C./matlab_normalization;
end

%% Pearson correlation [-1,1]
if mat_norm~=1 % pearson normalization, max value=1, min value=-1
    for k=-l:1:l % more than the zero time lag is requested

        if k<0 % negative time lag
            A=a(abs(k)+1:end); B=b(1:end-abs(k));
            mA=mean(A); mB=mean(B);
            
            C(counter) = (sum((A-mA).*(B-mB)))/sqrt( sum((A-mA).^2)*sum((B-mB).^2) );
            counter=counter+1;
        end
        if k==0 % 0 time lag
            C(counter) = (sum((a-mean(a)).*(b-mean(b)))) / sqrt( sum((a-mean(a)).^2) * sum((b-mean(b)).^2) );
            counter=counter+1;
        end
        if k>0 % positive time lag
            A=a(1:end-k); B=b(k+1:end);
            mA=mean(A); mB=mean(B);
            
            C(counter) = (sum((A-mA).*(B-mB)))/sqrt( sum((A-mA).^2)*sum((B-mB).^2) );
            counter=counter+1;
        end
    end

end
%%
C=fliplr(C); % flip to agree with matlab's xcorr
end