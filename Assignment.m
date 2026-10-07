clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=71.25;          % Settings of group 5 (Groep005.txt)
cFuel='Gasoline';                                                           % Fuel of this group (other choices check Sp.Name)
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
% Make a choice for which type of solution method you want to use.



%% [2-3] Compressor :: placeholder so this file runs on its own, replace by the compressor part (must give T3 and P3)-
P3 = P3overP2*P2;                                                           % Given pressure ratio
T3 = interp1(sair_a,TR,s2thermal+Rg*log(P3/P2));                            % Isentropic: thermal entropy rises by Rg*ln(P3/P2)

for i=1:NSp                                                                 %Check
    hi3air(i) = HNasa(T3,SpS(i));
end

h3 = Yair*hi3air';



%% [3-4] Combustor :: composition before and after the combustor
% Needs T3 and P3 from the compressor. Fuel and air both enter at T3.
sPart = 'Combustor';
mairrate = AF*mfurate;                       % Air mass flow, AF = mair/mfuel [kg/s]
mtotrate = mairrate+mfurate;                 % Mass conservation: flow leaving [kg/s]
Yreac = (AF*Yair+Yfuel)/(AF+1);              % Unburnt mixture: AF kg air + 1 kg fuel
% Complete combustion: CxHy + nuO2 O2 -> x CO2 + y/2 H2O
ElF = SpS(1).Elcomp;                         % Atoms in the fuel, order [O H C N Ar]
nC = ElF(3);nH = ElF(2);                     % x and y: number of C and H atoms
nuO2 = nC+nH/4;                              % Moles O2 needed per mole of fuel
AFstoi = nuO2*Mi(2)/Mi(1)/Yair(2);           % Stoichiometric air-fuel ratio [kg/kg]
phi = AFstoi/AF;                             % Equivalence ratio (<1: excess air)
nfuel = Yreac(1)/Mi(1);                      % Moles of fuel per kg unburnt mixture
Yprod = Yreac;                               % Burnt mixture, N2 does not react
Yprod(1) = 0;                                % All fuel is burnt
Yprod(2) = Yreac(2)-nuO2*nfuel*Mi(2);        % O2 that is left over
Yprod(3) = Yreac(3)+nC*nfuel*Mi(3);          % CO2 that is formed
Yprod(4) = Yreac(4)+nH/2*nfuel*Mi(4);        % H2O that is formed
Rgreac = Runiv*sum(Yreac./Mi);               % Gas constant unburnt mixture [J/kg/K]
Rgprod = Runiv*sum(Yprod./Mi);               % Gas constant burnt mixture [J/kg/K]


%% [3-4] Combustor :: thermodynamic computations
P4 = P3;                                     % Combustion at constant pressure
v3 = 0;v4 = 0;                               % Velocities in the engine are neglected
for i=1:NSp
    hi3(i) = HNasa(T3,SpS(i));               % hi(T3), formation enthalpy included
end
hreac = Yreac*hi3';                          % Enthalpy of the unburnt mixture at T3
h4 = hreac+0.5*v3^2-0.5*v4^2;                % Energy balance: adiabatic, no work
hprod_a = Yprod*hia';                        % Enthalpy burnt mixture for range of T
T4 = interp1(hprod_a,TR,h4);                 % Interpolate h4 on hprod_a to get T4
for i=1:NSp
    hi4(i) = HNasa(T4,SpS(i));
end
h4check = Yprod*hi4';                        % Check: must give h4 back
% Print to screen
fprintf('\nStage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,3,4);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T3,T4);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P3/kPa,P4/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v3,v4);
fprintf('---  H      -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',hreac/kJ,h4check/kJ);
fprintf('---  Composition (table 2)  ---------\n');
fprintf('%8s| %9.2f  (phi = %5.3f)\n','AF',AF,phi);
for i=[1 2 5 3 4]
    fprintf('%8s| %9.5f %9.5f  [kg/kg]\n',SpS(i).Name,Yreac(i),Yprod(i));
end
fprintf('%8s| %9.2f %9.2f  [J/kg/K]\n','Rg',Rgreac,Rgprod);


%%[4-5] Turbine

sPart = 'Turbine';                                                          % Defines word turbine.                            

compressorPower = mairrate*(h3-h2);                                         %Calculates how much power the compressor requires from the turbine. The compressor acts only on the air mass flow.

h5 = h4-compressorPower/mtotrate;

T5 = interp1(hprod_a,TR,h5);                                                %Searches the combustion/product enthalpy curve and finds the temperature that corresponds to h5.

for i=1:NSp                                                                 %Starts a loop through all five species: Gasoline, O2, CO2, H2O, and N2.
    si4(i) = SNasa(T4,SpS(i));
    si5(i) = SNasa(T5,SpS(i));
end

s4thermal = Yprod*si4';                                                     %Combines the entropy values of all species according to the combustion/product mass fractions.

s5thermal = Yprod*si5';                                                     %Same thing at state 5

lnP5overP4 = (s5thermal-s4thermal)/Rgprod;                                  %Total entropy stays constant. That gives the pressure ratio from the change in thermal entropy.

P5 = P4*exp(lnP5overP4);                                                    %Removes natural logarithm and gives the actual pressure ratio.

turbinePower = mtotrate*(h4-h5);                                            %Recalculates the turbine power as a check. It should be identical to compressorPower.

S4 = s4thermal-Rgprod*log(P4/Pref);                                         %Calculates the full specific entropy at state 4, including the pressure contribution.

S5 = s5thermal-Rgprod*log(P5/Pref);                                         %Full entropy at state 5. Since the turbine is isentropic, it should give approx. S4=S5.

fprintf('\nStage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,4,5);     %Print results...
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T4,T5);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P4/kPa,P5/kPa);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h4/kJ,h5/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S4/kJ,S5/kJ);
fprintf('%8s| %9.2f  [kW]\n','Compressor power',compressorPower/kJ);
fprintf('%8s| %9.2f  [kW]\n','Turbine power',turbinePower/kJ);