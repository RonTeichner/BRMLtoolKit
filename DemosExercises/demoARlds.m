
function demoARlds
%DEMOARLDS demo of fitting an autoregressive model to data
import brml.*
close all;
% generate from an AR LDS
L=10; % number of AR coefficients
T=4000; % number of timepoints
sigma2V=0.001; % output noise
sigma2H=0.0001; % AR coefficient transition noise

% generate some training data:
%r=rand; r2=rand;
omega1=rand; omega2=rand; % => max possible omega is 1 [rad] <=> max possible freq is 1/(2*pi) [hz] => fs > 1/pi
fs = 0.5; % [hz]
tVec=(1:T)*fs;
for i=1:T
	%v(t)=sin(0.5*r*tVec(t))+sin(0.5*r2*tVec(t)) + 0.001*randn; 
    station1(i) = sin(omega1(i)*tVec(i)); 
    station2(i) = sin(omega2(i)*tVec(i));
    omega1(i+1)=omega1(i)+sigma2H*randn; omega2(i+1)=omega2(i)+sigma2H*randn;            
end
Tx = station1 + station2;
v = Tx + sigma2V*randn(1,numel(Tx)); 	
omega1Mod = omega1 - omega1(1); omega2Mod = omega2 - omega2(1);
subplot(4,1,1); hold all; plot(tVec,(omega1Mod(1:numel(tVec))./(2*pi))./1e-3); plot(tVec,(omega2Mod(1:numel(tVec))./(2*pi))./1e-3); xlabel('sec'); ylabel('mhz'); legend('f1','f2'); title('diff from carrier');
subplot(4,1,2); plot(tVec,station1); hold all; plot(tVec,station2); plot(tVec,v); xlabel('sec'); title('time series')
legend('s1','s2','v');

% Learn the AR coefficients:
[f g]=ARlds(v,L,sigma2V,sigma2H);

%subplot(4,1,3);imagesc(f); title('filtered AR coefficients')
subplot(4,1,3); hold all; 
for i=1:L
    plot(tVec,f(i,:));
end
xlabel('sec'); title('filtered AR coefficients')
%subplot(4,1,4);imagesc(g);title('smoothed AR coefficients')
subplot(4,1,4); hold all; 
for i=1:L
    plot(tVec,g(i,:));
end
xlabel('sec'); title('smoothed AR coefficients')

% reconstructions:
vFromAR = zeros(size(v)); vFromAR(1:L) = v(1:L);
send_V_every_N_samples = 5; 
vFromAR_mixed = zeros(size(v)); vFromAR_mixed(1:L) = v(1:L);
for i=L+1:(T-1)
    
    B = v(i-L:i-1);        
    vFromAR(i) = B*f(:,i+1);
    
    if mod(i,send_V_every_N_samples) == 0
        vFromAR_mixed(i) = v(i);
    else
        B = vFromAR_mixed(i-L:i-1);
        vFromAR_mixed(i) = B*f(:,i+1);
    end
    
end
figure;
hold all; plot(tVec,v); plot(tVec,vFromAR); 
plot(tVec,vFromAR_mixed);
title('reconstruction attempts');
xlabel('sec'); legend('v','vFromAR','sending AR and some real samples');