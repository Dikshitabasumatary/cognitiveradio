%close all
clc
clear all
%Initialising the count to be equal to 1
count =1;
%initialising outage probability
outage_prob_count=0;
%initialising time index
time_index=0;
SNR_Occupied=[];
SNR_Empty=[];
C_Occupied=[];
Pout_Occ=[];
Throughput_Occ=[];
AvgCapacity_Occ=[];
Pout_Empty=[];
Throughput_Empty=[];
AvgCapacity_Empty=[];
%%% Simulation Parameteres
dSD= 2; % Distance between SS to SD
dSP = 4; % Distance between SS to PD
PL= 3; % Path loss between any two node
iterations= 1000;
C_Empty=[];
TIME=[];
PdB = 0:1:30; % Battery Power at SS
%%% Define the channels
hSD = exprnd(1/dSD^PL, [1, iterations]); % Channel realization between SS to SD
gSP = exprnd(1/dSP^PL, [1, iterations]); % Channel realization between SS to SP
for z = 1:1
    for i = 1:iterations
    %Checking channel occupancy for x Iterations, 
    while count<10
       if count ==1
            N1 = 0 ;
            %%Generating an exponential random number between 1 and 10
            d = makedist('Exponential','mu',5);
            td = d.truncate(1,10)
            x= td.random(1,1);
            N2 = x;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Assigning value of N1 and N2 at the successive iterations
        elseif count >=1
            %%%%%Part to generate Poission random number between 0 and 10%%
            lambda = (0+10)/2; % To put the mean in the middle of the range
            d = makedist('poisson', 'lambda', lambda);
        % Truncated Poisson distribution
        %td = d.truncate(0, 10)
        x= floor(d.random(1,1)); % Generate 1 random numbers in the range [0 10]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N1 = N2 + x;
        %New N2 becomes N1 + a exponential random number between 1 and 10
        %%Generating an exponential random number
        d = makedist('Exponential','mu',5);
        %td = d.truncate(1,10)
        x=floor(d.random(1,1));
        % x being a exponential random number between 1 and 10
        N2 = N1 + x;
       end
    %Checking channel after every 4 sec, you can change it as per your convenience
    channel_check_time =3 ;
    %incrementing counter at this stage
    count = count +1;
    %If the last time interval is divisible by channel check time then H1 else channel H0 for empty
    if mod(round(N2), channel_check_time)==0;
        disp('Channel Occupied')
        IdB=3; % Interference constraint at SS
        I=10^(IdB/10);
        No=1;
        SerTime = exprnd(5, [1, iterations]); %Service time of primary network
        %%% Secondary Rate
        Rs_T=2
        for j=1:length(Rs_T)
            Rs=Rs_T(j);
            for i=1:length(PdB)
                In_Var =i;
                P=10^(PdB(i)/10);
                Ps= min(P, I./gSP);
                SNR = Ps.*hSD/No;
                SNR_Occupied=[SNR_Occupied,SNR];
                C = log2(1+SNR); % Allowable capacity by secondary network in underlay
                %C_Occupied = [C_Occupied,C];
                Pout(In_Var)= sum(C<Rs)/iterations; % Outage probablity of the secondary network in underlay                
                T1(In_Var) = Rs*(1-Pout(In_Var)); % Throughput of the secondary network in underlay
                AvgCapacity(In_Var) = sum(C)/iterations;
            end
        end
        Pout_Occ=[Pout_Occ, Pout];
        Throughput_Occ=[Throughput_Occ, T1];
        AvgCapacity_Occ=[AvgCapacity_Occ,AvgCapacity];
        
    else
        No=1;
        %%% Secondary Rate of interweave
        Rs_T=2;       
        for j=1:length(Rs_T)
            Rs=Rs_T(j);
            for i=1:length(PdB)
               In_Var =i;
                P=10^(PdB(i)/10);
                SNR = P.*hSD/No;
                SNR_Empty=[SNR_Empty,SNR];
                C = log2(1+SNR); % Allowable capacity by secondary network in interweave
                %C_Empty = [C_Empty,C];
                Pout(In_Var)= sum(C<Rs)/iterations; % Outage probablity of the secondary network in interweave               
                T2(In_Var) = Rs*(1-Pout(In_Var)); % Throughput of the secondary network in interweave
                AvgCapacity(In_Var) = sum(C)/iterations;
            end
        end
        Pout_Empty=[Pout_Empty, Pout];
        Throughput_Empty=[Throughput_Empty, T2];
        AvgCapacity_Empty=[AvgCapacity_Empty,AvgCapacity];
    end
end 
    end
end
%Graph Analysis
%1. Outage Probability

x1=1:30
figure;
plot(x1,Pout_Occ(1:30),'b');
xlabel('PdB');
ylabel('Pout Occupied');
title('Pout Occupied vs PdB');
figure;
plot(x1,Pout_Empty(1:30),'r');
xlabel('PdB');
ylabel('Pout Empty');
title('Pout Empty vs PdB');
Pout_hybrid = Pout_Occ(1:30)+ Pout_Empty(1:30)
figure;
plot(x1,Pout_hybrid)
xlabel('PdB');
ylabel('Pout');
title('Pout Hybrid vs PdB');

%2. Channel Capacity

x2=1:30
figure;
plot(x2,AvgCapacity(1:30),'b');
xlabel('PdB');
ylabel('AvgCapacity when Channel is Occupied');
title('AvgCapacity Occupied vs PdB');
figure;
plot(x2,AvgCapacity_Empty(1:30),'r');
xlabel('PdB');
ylabel('AvgCapacity Empty');
title('AvgCapacity Empty vs PdB');
AvgCapacity_hybrid = AvgCapacity_Occ(1:30)+ AvgCapacity_Empty(1:30)
figure;
plot(x2,AvgCapacity_hybrid)
xlabel('PdB');
ylabel('AvgCapacity ');
title('AvgCapacity Hybrid vs PdB');



%3. Throughput
x=1:30
figure;
plot(x,Throughput_Occ(1:30),'b');
xlabel('PdB');
ylabel('Throughput Occupied');
title('Throughput Occupied vs PdB');
figure;
plot(x,Throughput_Empty(1:30),'r');
xlabel('PdB');
ylabel('Throughput Empty');
title('Throughput Empty vs PdB');
throughput_hybrid = Throughput_Occ(1:30)+ Throughput_Empty(1:30)
figure;
plot(x,throughput_hybrid)
xlabel('PdB');
ylabel('Throughput ');
title('Throughput Hybrid vs PdB');

% end
% % plot(T1,'y')
% % hold on;
% % plot(2*x^2* T1,'r')
% %Pout1=[Pout_Occ, Pout_Empt,];
% %Throughput1=[Throughput_Occ, Throughput_Empty];
% AvgCapacity1=[AvgCapacity_Occ,AvgCapacity_Empty];
% 
% %Graph Analysis
% figure
% plot(PdB, Pout, 'r');
% title('Outage Probability vs PdB of Interweave');
% xlabel('Pdb');
% ylabel('Outage Probability');
% figure;
% plot(PdB, T2,'b'); 
% title('ThroughPut vs PdB of Interweave');
% xlabel('Pdb');
% ylabel('ThroughPut');
% %hold on
% figure;
% plot(PdB,AvgCapacity,'g');
% title('Average Capacity vs PdB of Interweave');
% xlabel('Pdb');
% ylabel('Average Capacity');
