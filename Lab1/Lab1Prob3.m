clear all;
clc;

S0 = 100;
K=105;
T = 5;
r=0.05;
sigma = 0.3;

M=20;
dt = T/M;
    
    u = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    d = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
    p = (exp(r*dt)-d)/(u-d);
    
    if d<exp(r*dt) && exp(r*dt)<u
    %    disp('There is no arbitrage possible. Proceeding to calculate option prices');
    else
    %    disp('There is an arbitrage opportunity possible. The program will terminate');
        return;
    end
    
    [callprices, putprices]  =  tab(M,S0,K,r,u,d,p,dt);
    
calltab=table(callprices(:,1),callprices(:,3),callprices(:,5),callprices(:,7),callprices(:,13),callprices(:,19));
puttab=table(putprices(:,1),putprices(:,3),putprices(:,5),putprices(:,7),putprices(:,13),putprices(:,19));    
calltab.Properties.VariableNames = {'T0','T05','T10','T15','T30','T45'};
puttab.Properties.VariableNames = {'T0','T05','T10','T15','T30','T45'};

disp(calltab);
disp(puttab);

function [callprice, putprice]  =  tab(M,s0,K,r,u,d,p,dt)

    callprice = zeros(M+1,M+1);
    putprice = zeros(M+1,M+1);
    
    for i=1:M+1
        sn=d^(i-1)*u^(M-i+1)*s0;
        callprice(i,M+1) = max(0,sn-K);
        putprice(i,M+1) = max(0,K-sn);
    end
    
    for j = M:-1:1
        for  i =1:j
            callprice(i,j) = exp(-r*dt)*(p*callprice(i,j+1)+(1-p)*callprice(i+1,j+1));
            putprice(i,j) = exp(-r*dt)*(p*putprice(i,j+1)+(1-p)*putprice(i+1,j+1));
        end
    end

end