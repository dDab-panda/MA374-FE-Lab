clear all;
clc;


S0 = 100;
K=105;
T = 5;
r=0.05;
sigma = 0.3;

M1s = 1:1:100;
calloptions1 = zeros(1,length(M1s));
putoptions1 = zeros(1,length(M1s));

for i=1:length(M1s)
   
    dt = T/M1s(i);
    
    u = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    d = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    p = (exp(r*dt)-d)/(u-d);
    
    if d<exp(r*dt) && exp(r*dt)<u
    %    disp('There is no arbitrage possible. Proceeding to calculate option prices');
    else
    %    disp('There is an arbitrage opportunity possible. The program will terminate');
        return;
    end
    
    [calloptions1(i),putoptions1(i)] = cal(M1s(i),S0,K,r,u,d,p,dt);
end

figure(1);
plot(M1s,calloptions1,'r');
title('Call Options vs No. of Steps in binomial model(step size=1)');
xlabel('Number of Subintervals');
ylabel('Price of Call option');

figure(2);
plot(M1s,putoptions1,'b');
title('Put Options vs No. of Steps in binomial model(step size=1)');
xlabel('Number of Subintervals');
ylabel('Price of Put option');

M5s = 1:5:100;
calloptions2 = zeros(1,length(M5s));
putoptions2 = zeros(1,length(M5s));
for i=1:length(M5s)
   
    dt = T/M5s(i);
    
    u = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    d = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    p = (exp(r*dt)-d)/(u-d);

    [calloptions2(i),putoptions2(i)] = cal(M5s(i),S0,K,r,u,d,p,dt);
end

figure(3);
plot(M5s,calloptions2,'r');
title('Call Options vs No. of Steps in binomial model(step size=5)');
xlabel('Number of Subintervals');
ylabel('Price of Call option');

figure(4);
plot(M5s,putoptions2,'b');
title('Put Options vs No. of Steps in binomial model(step size=5)');
xlabel('Number of Subintervals');
ylabel('Price of Put option');

function [calloption, putoption]  =  cal(M,s0,K,r,u,d,p,dt)

    callprice = zeros(1,M+1);
    putprice = zeros(1,M+1);
    
    for i=1:M+1
        sn=d^(i-1)*u^(M-i+1)*s0;
        callprice(i) = max(0,sn-K);
        putprice(i) = max(0,K-sn);
    end
    
    for j = M:-1:1
        for  i =1:j
            callprice(i) = exp(-r*dt)*(p*callprice(i)+(1-p)*callprice(i+1));
            putprice(i) = exp(-r*dt)*(p*putprice(i)+(1-p)*putprice(i+1));
        end
    end
    
    calloption = callprice(1);
    putoption = putprice(1);
end