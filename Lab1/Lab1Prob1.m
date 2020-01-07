clear all;
clc;

M = [1 5 10 20 50 100 200 400];

S0 = 100;
K=105;
T = 5;
r=0.05;
sigma = 0.3;

calloptions = zeros(1,length(M));
putoptions = zeros(1,length(M));

for i=1:length(M)
   
    dt = T/M(i);
    
    u = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    d = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    p = (exp(r*dt)-d)/(u-d);

    if d<exp(r*dt) && exp(r*dt)<u
    %    disp('There is no arbitrage possible. Proceeding to calculate option prices');
    else
    %    disp('There is an arbitrage opportunity possible. The program will terminate');
        return;
    end
    
    [calloptions(i),putoptions(i)] = cal(M(i),S0,K,r,u,d,p,dt);
    
end

T = table(M',calloptions',putoptions');
T.Properties.VariableNames = {'M' , 'CallOptions','PutOptions'};
disp(T);

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

