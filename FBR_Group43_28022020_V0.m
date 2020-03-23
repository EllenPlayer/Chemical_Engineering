function dydt= FBR_Group43_28022020_V0(l, y, E, Dp) %y is a state vector
%--Defining elements in the matrix--
T= y(1);       % Temperature K
P= y(2);       % Pressure Bar
ni= y(3:7);    % Molar Flow rate kmol/s CO H M Me W

%--Calculating nt yi Pi--
nt=sum(ni);             % Total Molar Flow Rate kmol/s
yi=ni/nt;               % Composition of component i at a point
Pi=yi*P;                % Partial Pressure of component i(1:5) CO H M Me W
%--Data Given ---Constant mass flow rate in the reactor--
K=10^(6032/T - 10.531); % Rate constant from definition bar-2
% E=0.4;                  % Bed Voidage      (Commented out as defined when calling the function due to flexibility in changing the variable)
% Dp= 7.04* 10^(-3);      % Equivalent Diameter for assuming catalyist is a sphere m  
R= 0.08314;             % Gas Constant kJ/kmol
H1=-35490;              % Standard Heat of reaction 298K kJ/kmol
H2=-82518;              % Standard Heat of reaction 298K kJ/kmol
ParticleDencity=1400;   % kg/m3
Area=0.23;              % m2

%--Defining constants for reaction valid for range 550-750K--
A= -1103.6 -308.4 *(T/100) +171.6*(T/100)^2 -14.64*(T/100)^3;
B=92.65-17.95*(T/100)-0.5265*(T/100)^2+0.1745*(T/100)^3;
C=34.295 -5.444*(T/100)-0.6771*(T/100)^2 +0.1088*(T/100)^3;
D=483.5-28.21*(T/100) -22.8*(T/100)^2 +2.438*(T/100)^3;
%Reaction 1 CO + 2H2 -> CH3OH
R1= ((Pi(1)*(Pi(2))^2 -Pi(3)/K)/((A+B*Pi(1)+C*Pi(2)+D*Pi(3))^3))/3600; %kmolCO/kgcat.s
% Side Reaction 2 CO + 3H2 -> CH4 + H2O
R2 = (0.5*10^15 *exp(-30000/T - 1.6)*Pi(1)*(Pi(2)).^0.1)/3600; %kmolCO/kgcat.s
Reaction= [R1 ; R2];

%Stochiometric coefficents matrix CO H Me M W for R1 and R2
Stoc=[-1 -1; -2 -3; 1 0; 0 1; 0 1  ];
TransStoc1=Stoc(:,1)';%Transposing matrix for later multiplication
TransStoc2=Stoc(:,2)';%Transposing matrix for later multiplication

% --Defining heat capacity cp CO,H,M,Me,W--
CapasityPar=[ 42.2130 0.2505*10^(-2) 0.8055*10^(-5) -3.330*10^(-9) ;43.6305 -0.2865*10^(-2) 0.600*10^(-5) -1.3050*10^(-9) ; 28.5555 13.7190*10^(-2) -1.8255*10^(-5) -12.0495*10^(-9) ;29.8110 7.5315*10^(-2) 1.9020*10^(-5) -16.5000*10^(-9); 48.3255 0.2880*10^(-2) 1.5825*10^(-5) -5.3895*10^(-9)];
CpPure=CapasityPar*[1;T;T^2;T^3]; %Matrics of pure components Cpi KJ/KmolK
Cp = yi.' * CpPure; % Heat Capacity under specified conditions KJ/KmolK
cpPureint=CapasityPar *[T-298; (T^(2)-298^2)/2;(T^(3)-298^3)/3;(T^(4)-298^4)/4]; % intergrated Pure Cp values
CPMol=sum(ni*Cp); 
deltaHR1=H1+TransStoc1*cpPureint;  %Heats of reaction
deltaHR2=H2+TransStoc2*cpPureint;
%Molecular weights g/mol
MCO=28.01;
MH= 2.02;
MM= 32.04;
MMe= 16.04;
MW= 18.02;
ro = P*(yi(1)*MCO+yi(2)*MH + yi(3)*MM +yi(4)*MMe+yi(5)*MW)/(R*T); % Gas Density Kg/m3
vT = (ni(2)*MH + ni(3)*MM +ni(4)*MMe+ni(5)*MW+ni(1)*MCO)/ro;      %Volumetric Flow rate m3/s

%--Reactor Design Derivation Equations--
dTdl=ParticleDencity *(1-E)*Area *(-R1*(deltaHR1)-R2*(deltaHR2))/CPMol;
dPdl= (-1.75*(vT.^(2)) *ro *(1-E)/ ((Dp*Area^2 *(E)^3))*10^(-5));
dndl=Stoc * Reaction * Area * ParticleDencity * (1-E);
dydt=zeros(7,1) ; %Defining a vector for the elements
dydt=[dTdl ; dPdl; dndl]; %The state vetor elements returned as a fn of l
end