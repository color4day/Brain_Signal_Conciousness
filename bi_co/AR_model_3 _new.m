%% Autoregressive Modelling of Brain Signals
% Author: Kurt Butler, Yuanqing Song
% Date: Jan 5, 2021
% Rev: 3
% Description:
% Simulates the behavior of brain signals in the thalamus and PFC as
% modelled by simple AR and linear relationships. We model two thalamic
% channels and two cortex channels in order to understand how bipolar
% processing affects causal inference. In our case, we consider
% cross-correlation analysis as one measure of directionality. 

% General parameters
T0 = 10450;  % No. of samples in time (we pretend that each sample is 1 ms)
burn_in = 10; % No. of samples to be discarded at the beginning of a simulation

%% Parameter Initialization
T = T0 + burn_in; % Total no. of time points to be simulated

% AR coefficients
A = 0.97;
beta = [0.96 0.96 0.98 0.979];

% Scaling coefficients
a = 0*[1 1 1 1];% Common signal multiplier
b =  [1 1 1 1]; % Local signal multiplier
c = [NaN NaN 1 0; ... % Thalamus to PFC talk
     NaN NaN 1 0]; 

% Delays
tau = [0 0 0 0];
d = [0 0 10 10];

% Observation Noise Variance
sigw = [1 1 1 1];

% Covariance of local signal driving noise in the thalamus
thalamus_cov = 0.1;

%% Simulation
% Initialize 
s = zeros(T,4);
x = zeros(T,1);
y = zeros(T,4);
z = zeros(T,2);

% Local signal driving noise      thalamic noises are correlated     PFC noises are not
u = cat(2,  mvnrnd([0,0],toeplitz([1, thalamus_cov]),T),  randn(T,2) );  
% Common signal driving noise
v = normrnd(T,1);
% Observation noise
w = mvnrnd([0,0,0,0], diag(sigw),T);

% Signals
for t = 1:T
    for k = 1:4
        %% Neurons/local signals (S-Layer)
        % Driving noise u(k) 
        s(t,k) = u(t,k);
        
        % Autoregressive component (local effects)
        if t >= 2
            s(t,k) = s(t,k) + beta(k)*s(t-1,k);
        end
        
        %% Observed signals (Y-layer)
        % Local signals s(t) + Observation noise w(t) 
        y(t,k) = b(k)*s(t,k) + w(t,k)  ;
        
        if t > tau(k) % Common noise  x(t)  (common effects)
           y(t,k) = y(t,k) + a(k)*x(t); 
        end
        
        % Thalamic to PFC connections (topological effects)
        if k >= 3
            if t > d(1) 
                y(t,k) = y(t,k) + c(1,k)*s(t-d(1),1);
            end
            if t > d(2)
                y(t,k) = y(t,k) + c(2,k)*s(t-d(2),2);
            end
        end
    end
end

% Bipolar Processing  (Z-layer)
z(:,1) = y(:,2) - y(:,1);
z(:,2) = y(:,4) - y(:,3);

%% Plotting
time = (1:T)';
figure(1);
subplot(4,1,1);
plot(time,x,'r');
grid on, xlabel('Time (ms)'), title('X-process');
subplot(4,2,3);
plot(time,s(:,1),time,s(:,2));
grid on, title('Thalamic S-processes');
subplot(4,2,4);
plot(time,s(:,3),time,s(:,4));
grid on, title('PFC S-processes');
subplot(4,2,5);
plot(time,y(:,1),time,y(:,2));
grid on, title('Thalamic Y-processes');
subplot(4,2,6);
plot(time,y(:,3),time,y(:,4));
grid on, title('PFC Y-processes');
xlabel('Time');
subplot(4,2,7);
plot(time,z(:,1));
grid on, title('Thalamic Z-process');
subplot(4,2,8);
plot(time,z(:,2));
grid on, title('PFC Z-process');



%% Lagging Plots
LAGS = -50:50;
XCORR = zeros(size(LAGS,1),4);
for iter = 1:numel(LAGS)
    lag = LAGS(iter);
    id1 = max(1,1+lag):min(T,T+lag);
    id2 = max(1,1-lag):min(T,T-lag);
    XCORR(iter,1) = corr( z(id1,1), z(id2,2) );
end

% Plot it
figure(100);
plot(LAGS,XCORR(:,1))
grid on, title('Xcorr r_{z_1 z_2}(\tau)');
ylim([-1 1]);
