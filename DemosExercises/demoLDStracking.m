function demoLDStracking
%DEMOLDSTRACKING demo of tracking using a LDS Newtonian system
import brml.*
close all;
figure
T = 400; % number of timesteps
Delta = 0.1; % discretisation of time
% Hidden Variables are (in order) : xp,x,yp,y,fxom,fyom
% xp : x speed
% x : x position
% yp : y speed
% y : y position
% fxom : x acceleration
% fyom : y acceleration

% Newtonian Dynamics for Transitions:
A(1,:) = [1,0,0,0,Delta,0];
A(2,:) = [Delta,1,0,0,0,0];
A(3,:) = [0,0,1,0,0,Delta];
A(4,:) = [0,0,Delta,1,0,0];
A(5,:) = [0,0,0,0,1,0];
A(6,:) = [0,0,0,0,0,1];

% Observations are positions x and y
B(1,:) = [0,1,0,0,0,0];
B(2,:) = [0,0,0,1,0,0];

% Generate some Data:
h(2,1)=rand; h(4,1)=rand; % initial x and y position
h(1,1)=15*rand; h(3,1)=15*rand; % initial x and y velocity
h(5,1)=rand; h(6,1)=-rand; % initial x and y accelerations
sigV = 50;  % emission noise (standard deviation)-- same for both x and y
sigH = 0.00001; % small transition noise

v(:,1)= B*h(:,1)+sigV*randn(2,1); 
for t=2:T
  %h(:,t)=A*h(:,t-1)+sigH*randn(6,1); % Noisy Newtonian Dynamics  
  h(:,t)=A*h(:,t-1)+sigH*[zeros(4,1) ; randn(2,1)]; % Noisy Newtonian Dynamics  
  v(:,t)= B*h(:,t)+sigV*randn(2,1); % Noisy observation
end

% Setup the LDS :
stdFactorH = 1; stdFactorV = 1;
CovH=(stdFactorH*sigH)^2*eye(6); CovV=(stdFactorV*sigV)^2*eye(2); % transition and emission noise
CovP=1*eye(6); meanP=zeros(6,1); % vague prior
meanH=zeros(6,1);
meanV=zeros(2,1);
%Trajectory Estimates:
[dum1,dum2,mean_post,cov_post,dum3]=LDSsmooth(v,A,B,CovH,CovV,CovP,meanP,meanH,meanV,'cell');
hold on
for t=1:T;
	plot(v(1,t),v(2,t),'go','markersize',5);
    d = 10;
	if mod(t,d)==1 % only plot every 10th timestep (otherwise too messy)
		mh = mean_post{t};
		plot(h(2,t),h(4,t),'rx','markersize',10); plot(mh(2,:),mh(4,:),'b+')
        
        if ( (floor(t-d/2) >= 1) && (ceil(t+d/2) <= T))
            plot(mean(v(1,floor(t-d/2):ceil(t+d/2))) , mean(v(2,floor(t-d/2):ceil(t+d/2))) ,'m.','markersize',10);
            plot(median(v(1,floor(t-d/2):ceil(t+d/2))) , median(v(2,floor(t-d/2):ceil(t+d/2))) ,'k.','markersize',10);
        end
	end
end
xlabel('x','fontsize',10); ylabel('y','fontsize',10); set(gca,'box','on');

figure;
subplot(2,1,1);
hold on
for t=1:T
    x = h(2,t); 
	plot(t,v(1,t) - x,'go','markersize',5);
    d = 10;
	if mod(t,d)==1 % only plot every 10th timestep (otherwise too messy)
		mh = mean_post{t};
		plot(t,h(2,t) - x,'rx','markersize',10); plot(t,mh(2,:) - x,'b+')
        
        if ( (floor(t-d/2) >= 1) && (ceil(t+d/2) <= T))
            plot(t,mean(v(1,floor(t-d/2):ceil(t+d/2))) - x ,'m.','markersize',10);
            plot(t,median(v(1,floor(t-d/2):ceil(t+d/2))) - x, 'k.','markersize',10);
        end
	end
end
xlabel('t','fontsize',10); ylabel('x','fontsize',10); set(gca,'box','on');
subplot(2,1,2);
hold on
for t=1:T
    y = h(4,t); 
	plot(t,v(2,t) - y,'go','markersize',5);
    d = 10;
	if mod(t,d)==1 % only plot every 10th timestep (otherwise too messy)
		mh = mean_post{t};
		plot(t,h(4,t) - y,'rx','markersize',10); plot(t,mh(4,:) - y,'b+')
        
        if ( (floor(t-d/2) >= 1) && (ceil(t+d/2) <= T))
            plot(t,mean(v(2,floor(t-d/2):ceil(t+d/2))) - y ,'m.','markersize',10);
            plot(t,median(v(2,floor(t-d/2):ceil(t+d/2))) - y, 'k.','markersize',10);
        end
	end
end
xlabel('t','fontsize',10); ylabel('y','fontsize',10); set(gca,'box','on');