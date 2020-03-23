%% Reactor Design and Control Project - Group 43
% Ellen Player, Xiaozhou Feng, Mia McLachlan
%---------------------------------------------------------------------
%% Clearing previous imputs to avoid interference from historic data
clear;clc;
%% Obtaining the required graphs

E=0.4;                  % Bed Voidage 
Dp= 7.04* 10^(-3);      % Equivalent Diameter for assuming catalyist is a sphere m
F0=1; % Feed Rate kmol/s
yo=[600 ; 450; F0*0.53;F0*0.43;0;F0*0.02;F0*0.02 ]; % Matrix of inital Conditions T(K) P(Bar) Component Flow Rates (Kmol/s) CO H M Me W  
lspan=[0 120]; %Investigsting FBR length l for a resonable range of l 0-120m
[l,y]= ode45(@(l, y) FBR_Group43_28022020_V0(l, y, E, Dp) ,lspan,yo); % Solving the simultaneous equation relative to l and inserting into a matrix length and variables. Solving for specified E and Dp

% tiledlayout(3,1) %To test l T P F relationships
% nexttile
% hold on
figure(1) % Comment out this if testing the T P F relationships in tiled layout
hold on
plot(l,y(:,3:7))
xlabel('length of reactor m')
ylabel('Molar Flow Rate kmol')
plot([0 max(l)],[0.077 0.077]) % This line plotted on the graph shows the specified production capacity of the plant of 0.077kmol/s in the plant.
legend({'CO','H2','M','Me','W'})
hold off
% nexttile
% plot(l,y(:,1))
% xlabel('length of reactor m')
% ylabel('Temperature within the FBR K')
% nexttile
% plot(l,y(:,2))
% xlabel('length of reactor m')
% ylabel('Presure in the FBR Bar')
% hold off
figure(2)  % Plotting the extent of reaction 1 relationships with other variables in the reactor T,P,Catalyist Weight,Extent of reaction 2
tiledlayout(4,1) % Keeps graphs together in a column
nexttile
hold on
XI1=y(:,5)-yo(5); %Extent of reaction for reaction 1 kmol/s
XI2=y(:,6)-yo(6); %Extent of reaction for reaction 2 kmol/s
Area=0.23;              % m2
ParticleDencity=1400;   % kg/m3
CatalyistWeight=Area*l*(1-E)*ParticleDencity; %Kg
plot(XI1,CatalyistWeight,'-');
axis tight
grid on
title ('Relationship of Catalyist weight and Extent of reaction 1 with variable conditions')
xlabel ('Extent of Reaction 1 (kmol/s)')
ylabel ('Catalyist weight required (kg)')
nexttile
plot(XI1,XI2,'-');
axis tight
grid on
title ('Observing the reaction progression play off in the extents of the two reactions')
xlabel ('Extent of Reaction 1 (kmol/s) ')
ylabel ('Extent of Reaction 2 (kmol/s)')
nexttile
plot(XI1,y(:,1),'-');
axis tight
grid on
title ('Relationship of the Temperature as a result of extent')
xlabel ('Extent of Reaction 1 (kmol/s)')
ylabel ('Temperature (K)')
nexttile
plot(XI1,y(:,2),'-');
axis tight
grid on
title ('Relationship of the Pressure as a result of extent')
xlabel ('Extent of Reaction 1 (kmol/s)')
ylabel ('Pressure (bar)')
hold off
%% Obtaining additional graph of Reaction Rate vs Extent for the initial reaction

ni= y(:,3:7);   %Retreiving the molar flow rate of a given component at any given length in the reactor. kmol/s
nt=sum(ni,2);   %Calculating the total molar flow rate at a given point. kmol/s
yi=zeros(49,5); %Molar composition
pi=zeros(49,5); %Partial Pressure
for i= 1:5  % For loop calculating y and p for each component i
    yi(:, i)=ni(:,i)./nt;
    pi(:,i)=yi(:,i).*y(:,2);
end
T= y(:,1);  % Matrix of temperatures in the reactor for a given length
K=10.^(6032./T - 10.531); %Equilibrium Constant Bar^(-2)
%--Defining constants for reaction valid for range 550-750K--
A= -1103.6 -308.4*(T./100) +171.6*(T./100).^2 -14.64*(T./100).^3;
B=92.65-17.95*(T./100)-0.5265*(T./100).^2+0.1745*(T./100).^3;
C=34.295 -5.444*(T./100)-0.6771*(T./100).^2 +0.1088*(T./100).^3;
D=483.5-28.21*(T./100) -22.8*(T./100).^2 +2.438*(T./100).^3;
%Reaction 1 CO + 2H2 -> CH3OH
R1= ((pi(:,1).*(pi(:,2)).^2 -pi(:,3)./K)./((A+B.*pi(:,1)+C.*pi(:,2)+D.*pi(:,3)).^3))./3600;  % Rate of reaction 1 equation
figure(3)
plot(XI1, R1)
axis tight
grid on
title ('Rate vs Extent')
xlabel ('Extent of Reaction 1 (kmol/s) ')
ylabel ('Rate of Reaction 1(kmol/s)')
%% Analyisis of varying F0 for a fixed reactor length of chosen value- We assesed 120m 

%The code allows you to imput your own value of length dependent on the
%fixed sizes avaliable to be brought the number will need to  be inputed
%into the command window
prompt = 'What is the reactor length (m)? ';
Reactorlength = input(prompt);   %Used as 120 in our analysis

CatalyistWeight1=Area*Reactorlength*(1-E)*ParticleDencity; %kg
Q=zeros(5,4); %To increase inputs in i change the 5 value to the number of inputs added
c = 0; % initial value of c to run the first iteration
for i=[ 1 1.1 1.2 1.3 1.4 ] % Picking a range of initial flowrates- this range is chosen as outside of this pressure and temperature effects exceed the reactors capability
    c = c + 1;
    yo=[600 ; 450; i*0.53;i*0.43;0;i*0.02;i*0.02 ]; % Matrix of inital Conditions T(K) P(Bar) Component Flow Rates (kmol/s) CO H M Me W  
    lspan=[0 Reactorlength]; %Investigsting FBR length l for a resonable range of l 0-120m
    [l,y]= ode45(@(l,y) FBR_Group43_28022020_V0(l,y,E,Dp) ,lspan,yo);
    Q(c, :) = [y(end,5)  y(end,1)  y(end,5)/CatalyistWeight1 y(end,2)];
end 
figure(4) %Plotting a 3 axis graph
hold on 
x=Q(:,1);
y=Q(:,2);
z=Q(:,3);
r=Q(:,4);
yyaxis left
plot(x,y)
xlabel('Methanol Flow rate (kmol/s)')
ylabel('Temperature in reactor (T)')
yyaxis right
plot(x,z)
ylabel('Efficiency of catalyist use (kmol/kgs)')
axis tight
grid on
hold off
figure(5) %Observing the pressure drop for further analysis
plot(x,r)
xlabel('Methanol Flow rate (kmol/s)')
ylabel('Pressure (Bar)')

%% Abbreviations/Nomenclature
% CO= Carbon Monoxide
% H = Hydrogen
% M = Methanol
% Me= Methane
% W = Water