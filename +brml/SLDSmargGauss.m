function [meanGauss, covGauss]=SLDSmargGauss(weights,ps,means,covs)
%SLDSMARGGAUSS compute the single Gaussian from a weighted SLDS mixture
% [meanGauss, covGauss]=SLDSmargGauss(weights,ps,means,covs)
[S T]=size(ps); % S - no. of possible states; T - no. of time steps
J=size(means,2); % J - no. of gaussian components per state at each time step
for t=1:T;	
    ind=0;
    for s=1:S
        for j=1:J
            ind=ind+1; coeff(ind)=weights(j,s,t)*ps(s,t);
            mn(:,ind) = means(:,j,s,t); cv(:,:,ind)=covs(:,:,j,s,t);
        end
        % for a single time-step we take all the gaussians (all states and
        % all components per state) and collapse it to a single gaussian. so
        % we remain with a single gaussian per time-step:
        [dum,meanGauss(:,t),covGauss(:,:,t)]=brml.mix2mix(coeff,mn,cv,1); % find the mean of the Gaussian mixture
    end
end
