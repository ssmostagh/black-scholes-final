k = 25;
sigma = 0.18;
r = 0.0;
% T = [1:7];
s0 = 20;
V = zeros(7001,1);
i = 1;
mu = 0.15;

randn('state',100); % set the state of randn
T = 1; N = 7001; dt = 0.001;
dW = zeros(1,N); % preallocate arrays ...
W = zeros(1,N); % for efficiency
dW(1) = sqrt(dt)*randn; % first approximation outside the loop ...
W(1) = dW(1); % since W(0) = 0 is not allowed

for j = 2:N
dW(j) = sqrt(dt)*randn; % general increment
W(j) = W(j-1) + dW(j);
end

for T = 0:0.001:7
    s = s0*exp(sigma*W(i) + mu*T);
    x = ((log(s/k) + (r + (sigma^2/2))*T)/(sigma*sqrt(T)));
    z = ((log(s/k) + (r - (sigma^2/2))*T)/(sigma*sqrt(T)));
    PHI = @(j) (1/sqrt(2*pi))*integral(@(y) exp(-y.^2/2),-Inf,j);
    V(i) = s*PHI(x) - k*exp(-r*T)*PHI(z);
    i = i+1;
end

plot(V,'y')
title('Black-Scholes with Brownian Motion')
xlabel('time')
ylabel('Portfolio Value')