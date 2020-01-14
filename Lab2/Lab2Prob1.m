clear all;
clc;

M = 100;
S0 = 100;
K=100;
T = 1;
r=0.08;
sigma = 0.2;

dt = T/M;
u1 = exp(sigma*sqrt(dt));
d1 = exp(-sigma*sqrt(dt));
p1= (exp(r*dt)-d1)/(u1-d1);
u2 = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
d2 = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
p2 = (exp(r*dt)-d2)/(u2-d2);


% Question 1 Part (a)
% Varying S(0) and keeping other constant

S0 = 50:10:200;

callset1 = zeros(1,length(S0));
putset1 = zeros(1,length(S0));
callset2 = zeros(1,length(S0));
putset2 = zeros(1,length(S0));

for i=1:length(S0)
    [callset1(i),putset1(i)] = cal(M,S0(i),K,r,u1,d1,p1,dt);
    [callset2(i),putset2(i)] = cal(M,S0(i),K,r,u2,d2,p2,dt);
end

figure ;
subplot(2,2,1);
plot(S0,callset1,'r');
title('Call Options vs Price of Initial Stock with set1');
xlabel('Price of Initial Stock');
ylabel('Price of Call option');
subplot(2,2,2);
plot(S0,putset1,'g');
title('Put Options vs Price of Initial Stock with set1');
xlabel('Price of Initial Stock');
ylabel('Price of Put option');
subplot(2,2,3);
plot(S0,callset2,'r');
title('Call Options vs Price of Initial Stock with set2');
xlabel('Price of Initial Stock');
ylabel('Price of Call option');
subplot(2,2,4);
plot(S0,putset1,'g');
title('Put Options vs Price of Initial Stock with set2');
xlabel('Price of Initial Stock');
ylabel('Price of Put option');

% Question 1 Part (b)
% Varying K and keeping other constant


M = 100;
S0 = 100;
T = 1;
r=0.08;
sigma = 0.2;
K = 50:5:150;

callset1 = zeros(1,length(K));
putset1 = zeros(1,length(K));
callset2 = zeros(1,length(K));
putset2 = zeros(1,length(K));

for i=1:length(K)
    [callset1(i),putset1(i)] = cal(M,S0,K(i),r,u1,d1,p1,dt);
    [callset2(i),putset2(i)] = cal(M,S0,K(i),r,u2,d2,p2,dt);
end
figure ;
subplot(2,2,1);
plot(K,callset1,'r');
title('Call Options vs Different value of K with set1');
xlabel('Different value of K');
ylabel('Price of Call option');
subplot(2,2,2);
plot(K,putset1,'g');
title('Put Options vs Different value of K with set1');
xlabel('Different value of K');
ylabel('Price of Put option');
subplot(2,2,3);
plot(K,callset2,'r');
title('Call Options vs Different value of K with set2');
xlabel('Different value of K');
ylabel('Price of Call option');
subplot(2,2,4);
plot(K,putset1,'g');
title('Put Options vs Different value of K with set2');
xlabel('Different value of K');
ylabel('Price of Put option');


% Question 1 Part (c)
% Varying r and keeping other constant


r=0.01:0.01:0.2;
M = 100;
S0 = 100;
K=100;
T = 1;
sigma = 0.2;

callset1 = zeros(1,length(r));
putset1 = zeros(1,length(r));
callset2 = zeros(1,length(r));
putset2 = zeros(1,length(r));

for i=1:length(r)
    u1 = exp(sigma*sqrt(dt));
    d1 = exp(-sigma*sqrt(dt));
    u2 = exp(sigma*sqrt(dt)+(r(i)-0.5*sigma*sigma)*dt);
    d2 = exp(-sigma*sqrt(dt)+(r(i)-0.5*sigma*sigma)*dt);
    p1= (exp(r(i)*dt)-d1)/(u1-d1);
    p2= (exp(r(i)*dt)-d2)/(u2-d2);
    
    [callset1(i),putset1(i)] = cal(M,S0,K,r(i),u1,d1,p1,dt);
    [callset2(i),putset2(i)] = cal(M,S0,K,r(i),u2,d2,p2,dt);
end
figure ;
subplot(2,2,1);
plot(r,callset1,'r');
title('Call Options vs Different value of r with set1');
xlabel('Different value of r');
ylabel('Price of Call option');
subplot(2,2,2);
plot(r,putset1,'g');
title('Put Options vs Different value of r with set1');
xlabel('Different value of r');
ylabel('Price of Put option');
subplot(2,2,3);
plot(r,callset2,'r');
title('Call Options vs Different value of r with set2');
xlabel('Different value of r');
ylabel('Price of Call option');
subplot(2,2,4);
plot(r,putset1,'g');
title('Put Options vs Different value of r with set2');
xlabel('Different value of r');
ylabel('Price of Put option');

% Question 1 Part (d)
% Varying sigma and keeping other constant


sigma=0.01:0.01:0.3;

M = 100;
S0 = 100;
K=100;
T = 1;
r = 0.08;

callset1 = zeros(1,length(sigma));
putset1 = zeros(1,length(sigma));
callset2 = zeros(1,length(sigma));
putset2 = zeros(1,length(sigma));

for i=1:length(sigma)
    u1 = exp(sigma(i)*sqrt(dt));
    d1 = exp(-sigma(i)*sqrt(dt));
    u2 = exp(sigma(i)*sqrt(dt)+(r-0.5*sigma(i)*sigma(i))*dt);
    d2 = exp(-sigma(i)*sqrt(dt)+(r-0.5*sigma(i)*sigma(i))*dt);
    p1= (exp(r*dt)-d1)/(u1-d1);
    p2= (exp(r*dt)-d2)/(u2-d2);
    
    [callset1(i),putset1(i)] = cal(M,S0,K,r,u1,d1,p1,dt);
    [callset2(i),putset2(i)] = cal(M,S0,K,r,u2,d2,p2,dt);
end
figure ;
subplot(2,2,1);
plot(sigma,callset1,'r');
title('Call Options vs Different value of sigma with set1');
xlabel('Different value of sigma');
ylabel('Price of Call option');
subplot(2,2,2);
plot(sigma,putset1,'g');
title('Put Options vs Different value of sigma with set1');
xlabel('Different value of sigma');
ylabel('Price of Put option');
subplot(2,2,3);
plot(sigma,callset2,'r');
title('Call Options vs Different value of sigma with set2');
xlabel('Different value of sigma');
ylabel('Price of Call option');
subplot(2,2,4);
plot(sigma,putset1,'g');
title('Put Options vs Different value of sigma with set2');
xlabel('Different value of sigma');
ylabel('Price of Put option');


% Question 1 Part (c)
% Varying r and keeping other constant



M = 50:1:150;
S0 = 100;
K= 95:5:105;
T = 1;
sigma = 0.2;
r = 0.08;


for j=1:length(K)
    
callset1 = zeros(1,length(M));
putset1 = zeros(1,length(M));
callset2 = zeros(1,length(M));
putset2 = zeros(1,length(M));

for i=1:length(M)
dt=T/M(i);
u1 = exp(sigma*sqrt(dt));
d1 = exp(-sigma*sqrt(dt));
p1= (exp(r*dt)-d1)/(u1-d1);
u2 = exp(sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
d2 = exp(-sigma*sqrt(dt)+(r-0.5*sigma*sigma)*dt);
p2 = (exp(r*dt)-d2)/(u2-d2);
    [callset1(i),putset1(i)] = cal(M(i),S0,K(j),r,u1,d1,p1,dt);
    [callset2(i),putset2(i)] = cal(M(i),S0,K(j),r,u2,d2,p2,dt);
end
figure ;
subplot(2,2,1);
plot(M,callset1,'r');
title('Call Options vs Different value of M with set1');
xlabel('Different value of M');
ylabel('Price of Call option');
subplot(2,2,2);
plot(M,putset1,'g');
title('Put Options vs Different value of M with set1');
xlabel('Different value of M');
ylabel('Price of Put option');
subplot(2,2,3);
plot(M,callset2,'r');
title('Call Options vs Different value of M with set2');
xlabel('Different value of M');
ylabel('Price of Call option');
subplot(2,2,4);
plot(M,putset1,'g');
title('Put Options vs Different value of M with set2');
xlabel('Different value of M');
ylabel('Price of Put option');

    
end


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


