function [inspk] = spike_sort(spikes,method,inputs)
scales = 5;
%inputs = 20;
nspk=size(spikes,1);
ls = size(spikes,2);

% CALCULATES FEATURES
switch method
    case 'wav'
        cc=zeros(nspk,ls);
        for i=1:nspk                                % Wavelet decomposition
            [c,l]=wavedec(spikes(i,:),scales,'haar');
            cc(i,1:ls)=c(1:ls);
        end
        for i=1:ls                                  % KS test for coefficient selection   
            thr_dist = std(cc(:,i)) * 3;
            thr_dist_min = mean(cc(:,i)) - thr_dist;
            thr_dist_max = mean(cc(:,i)) + thr_dist;
            aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i);
            
            if length(aux) > 10;
                [ksstat]=test_ks(aux);
                sd(i)=ksstat;
            else
                sd(i)=0;
            end
        end
        [max ind]=sort(sd);
        coeff(1:inputs)=ind(ls:-1:ls-inputs+1);
    case 'pca'
        [C,S,L] = princomp(spikes);
        cc = S;
        %inputs = 5; 
        coeff(1:inputs)=meshgrid(1:inputs,1);
end
inspk=zeros(nspk,inputs);
for i=1:nspk
    for j=1:inputs
        inspk(i,j)=cc(i,coeff(j));
    end
end
end

function [KSmax] = test_ks(x)
% 
% Calculates the CDF (expcdf)
[y_expcdf,x_expcdf]=cdfcalc(x);

%
% The theoretical CDF (theocdf) is assumed to be normal  
% with unknown mean and sigma

zScores  =  (x_expcdf - mean(x))./std(x);
theocdf  =  normcdf(zScores , 0 , 1);

%
% Compute the Maximum distance: max|S(x) - theocdf(x)|.
%

delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
deltacdf  =  abs([delta1 ; delta2]);

KSmax =  max(deltacdf);
end