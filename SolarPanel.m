% 09/24/2018
% Esteban Perez, Hamza Bakheet, and Ali Alaraini
% ME 499-008: Solar Panel Project

% This code will produce and export individual excel tables for all variables 
% that depend on days and hours, and it will produce a table for just the 
% variables that depend on days
% All results in the tables are in ST 
%
% This code will also produce different plots for each variable
% 
% This code will ask you for a month, a day, and an hour CT or ST, and it will ask
% you to choose a variable that you'd like to calculate. 
% 
% if you would like to see a different result, just use the command window
% by typing the variable after the code has been run.

%Variables
delta = 0;  %Solar declination
n = 0;   %day number
L = 0;   %latitude
beta = 0; %solar altitude angle
betaN = 0; %angle at solar noon
H = 0;      %hour anlge
Hsr = 0;    %sunrise hour angle
phis = 0;   %solar azimuth
phic = 0;   %collector azimuth angle
ST = 0;     %soalr time
CT = 0;     %civil clock time
E = 0;      %equation of time
sigma = 0;  %collector til angle
tetha = 0;  %incidence angle between sun and colector phace
hour = 0;
m = 0 ;     %air mas ratio
IB = 0;     %beam insolation at earth's surface
A = 0;      %apparent extraterrestrial solar insolation
k = 0;      %atmospheric optical depth
C = 0;      %sky diffuse factor
IBC = 0;    %beam insulation on collector
IDC = 0;    %diffuse insolation on collector
IRC = 0;    %reflected insulation on collector
rho = 0;    %ground reflectance
IC = 0;     %insolation on collector

clear
clc

%Location: Rappahannock River Parking Deck
%Latitude 38.834703,Longtitude -77.305833

format bank

L = 38.834703;
long = 77.305833;
phic = 0*ones(365,24); % phic = 0  %<--corrected
rho = .2*ones(365,24);

%Preallocating
H = zeros(1,24);
ST = H;
delta = zeros(1,365);
B = delta;
A = B;
C = A;
k = A;
E = A;
%sigma = A;
hsr =A;
hss =A;

%Daily Values
for n = 1:365
    delta(n) = 23.5*sind(360*(n-81)/365);
    B(n) = 360*(n-81)/364;
    A(n) = 1160 + 75*sind(360*(n-275)/365);
    C(n) = 0.095 + 0.04*sind(360*(n-100)/365);
    E(n) = 9.87.*sind(2.*B(n))-7.53.*cosd(B(n))-1.5.*sind(B(n));
    k(n) = 0.174 + 0.035*sind(360*(n-100)/365);
    %sunrise hour
    hsr(n) = 12-acosd(-tand(L)*tand(delta(n)))/15;
    %sunset hour
    hss(n) = acosd(-tand(L)*tand(delta(n)))/15+12;
    %sigma(n) = -1*(-L+delta(n));  %sigma perpendicular to sun
    
    
    N(n) = n;
end

%Solar Time
for h = 1:24
        if h >= 12
            H(h) = -15*(h-12);
        else
            H(h) = 15*(12-h);
        end  
        ST(h) = h;
end

de = transpose(delta);
k =  transpose(k);
%sigma = transpose(sigma);
sigma = 40;
A = transpose(A);
C = transpose(C);
E = transpose(E);
B = transpose(B);
beta = asind(cosd(L)*cosd(de)*cosd(H)+sind(L)*sind(de));
beta = max(beta,0);
phis1 = asind((cosd(de)*sind(H))./(cosd(beta)));
phis2 = 180 - asind((cosd(de)*sind(H))./(cosd(beta)));

xx = cosd(H);
xx2 = tand(de)/tand(L);

%solar azimuth angle 
for i = 1:365
    for j = 1:24
        if xx(1,j)>=xx2(i,1)
            phis(i,j)=phis1(i,j);
        else
            phis(i,j)=phis2(i,j);
        end
    end
end

m = (((708.*sind(beta)).^2)+1417).^.5-(708.*sind(beta));
%m= 1./sind(beta);
m =max(m,0);
tetha = acosd(cosd(beta).*cosd(phis-phic).*sind(sigma) + sind(beta).*cosd(sigma));
IB = A.*exp(-k.*m);
IBC = IB.*cosd(tetha);
IDC = IB.*C.*(.5+cosd(sigma)./2);
IRC = rho.*IB.*(sind(beta)+C.*ones(1,24)).*(.5-cosd(sigma)./2); %<--corrected
IC = IBC + IDC + IRC;

%Part 2 - Graphs
dd = (1:365)';
tf = { 'Day', 'delta', 'A', 'B', 'C', 'k', 'E'};
w = [ dd de A B C k E ];
w = [ tf ; num2cell(w)];
plot(N,B)
title('B vs Day (N)')
xlabel('Day')
ylabel('B')
figure
plot(N,A)
title('Apparent Extraterrestrial Flux(A) vs Day(N)')
xlabel('Day Number')
ylabel('Apparent Extraterrestrial Flux (W/m^2)')
figure
plot(N,C)
title('Sky Diffuse Factor(C) vs Day(N)')
xlabel('Day Number')
ylabel('Sky Diffuse Factor')
figure
plot(N,E)
title('Equation of Time(E) vs Day(N)')
xlabel('Day Number')
ylabel('Equation of Time (min)')
figure
plot(N,k)
title('Optical Depth(k) vs Day(N)')
xlabel('Day Number')
ylabel('Optical Depth')
figure
plot(N,de)
title('Solar Declination(de) vs Day(N)')
xlabel('Day Number')
ylabel('Solar Declination (Degrees)')

%Part 3 - Contour Plots
contourf(beta)
title('Day(N) vs Hour(hr) vs Solar Altitude Angle(beta)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(phis)
title('Day(N) vs Hour(hr) vs Solar Azimuth Angle(phis)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(m)
title('Day(N) vs Hour(hr) vs Air Mass Ratio(m)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(IBC)
title('Day(N) vs Hour(hr) vs Beam Insolation on Collector(IBC)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(IDC)
title('Day(N) vs Hour(hr) vs Diffuse Insolation on Collector(IDC)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(IRC)
title('Day(N) vs Hour(hr) vs Reflected Insolation on Collector(IRC)')
xlabel('Hour')
ylabel('Day Number')
colorbar
figure
contourf(IC)
title('Day(N) vs Hour(hr) vs Insolation on Collector(IC)')
xlabel('Hour')
ylabel('Day Number')
colorbar

%Table Labels 
HourNum = {'Days','1 am', '2 am','3 am', '4 am', '5 am', '6 am', '7 am', '8 am', '9 am', '10 am', '11 am', '12 pm','1 pm', '2 pm','3 pm', '4 pm', '5 pm', '6 pm', '7 pm', '8 pm', '9 pm', '10 pm', '11 pm', '12 am'};
beta2 = [dd beta];
beta2 = [ HourNum ; num2cell(beta2)];
m2 = [dd m];
tetha2 = [dd tetha];
IB2 = [dd IB];
IBC2 = [dd IBC];
IDC2 = [dd IDC];
IRC2 = [dd IRC];
IC2 = [dd IC];
m2 = [ HourNum ; num2cell(m2)];
tetha2 = [ HourNum ; num2cell(tetha2)];
IB2 = [ HourNum ; num2cell(IB2)];
IBC2 = [ HourNum ; num2cell(IBC2)];
IDC2 = [ HourNum ; num2cell(IDC2)];
IRC2 = [ HourNum ; num2cell(IRC2)];
IC2 = [ HourNum ; num2cell(IC2)];
phis2 = [dd phis];
phis2 = [ HourNum ; num2cell(phis2)];

%Excel Table Outputs for Part 2 ( DALY VALUES) AND Part 3 (HOURLY VALUES)
% if this does not work for any reason, call any of the variables in the
% command window like IC2, phis2..etc to see the tables in Matlab
xlswrite('Daily Values.xlsx',w)
xlswrite('Solar Altitude Angle.xlsx',beta2)
xlswrite('Solar Azimuth Angle.xlsx',phis2);
xlswrite('Air Mass Ratio.xlsx',m2);
xlswrite('Beam Insulation on Collector.xlsx',IBC2);
xlswrite('Diffuse Insulation on Collector.xlsx',IDC2);
xlswrite('Reflected Altitude Angle.xlsx',IRC2);
xlswrite('Insolation on Collector.xlsx',IC2);

% Part 4 of the project, this is ONLY calculating values for a SPECIFIC time 
% no tables are genarate here. 
mo = 'month: ';
d = 'day: ';

month = input(mo);
day = input(d);
 
if month == 1
    n = day;
elseif month == 2
    n = day + 31;
elseif month == 3
    n = day + 59;
elseif month == 4
    n = day + 90;
elseif month == 5
    n = day + 120;
elseif month == 6
    n = day + 151;
elseif month == 7
    n = day + 181;
elseif month == 8
    n = day + 212;
elseif month == 9
    n = day + 243;
elseif month == 10
    n = day + 273;
elseif month == 11
    n = day + 304;
elseif month == 12
    n = day + 334;
end

de = delta(n);
B = B(n);
A = A(n);
C = C(n); 
k = k(n);

%calculating the exact Sunrise and Sunset time
hr2 = hsr(n);
hs = hss(n);
Q = 3.467/(cosd(L)*cosd(de)*sind(12*15-hr2*15));
E = E(n);
hr = hr2 + (( - 4*(75-long)) - E)/60-Q/60;
hs = hs + (( - 4*(75-long)) - E)/60+Q/60;
mhsr = abs(floor(hr)-hr)*60;
mhss = abs(floor(hs)-hs)*60;

fprintf('Sunrise at (CT): %d:%f am \n', floor(hr),mhsr);
fprintf('Sunset at (CT): %d:%f  pm \n', floor(hs)-12,mhss);
fprintf('day selcted: %d\n', n);

%Time
time = 'select a time with minutes as hour fractions(e.g 1:30 = 1.5): ';
t = input(time);
pam = 'pm[1] or am[0](select 0 for noon):  ';
pa = input(pam);
if pa == 1
    t = (t+12);
else
    
end

za = ' Is your input CT or ST? CT[0] ST[1]: ';
za =  input(za);
if za == 0
     t2 = t;
     t = t + (( 4*(75-long)) + E)/60;
    if t >= 25
        t = t-24;
        n = n+1;
    end
        fprintf('Solar Time Slected: %f hours\n', t)
        fprintf('Civil Time Slected: %f hours\n', t2)
elseif za == 1
       t2 = t - (( 4*(75-long)) - E)/60;
       fprintf('Solar Time Slected: %f hours\n', t)
       fprintf('Civil Time Slected: %f hours\n', t2)
end
 
%if the solar time is 25, it means that is the next day at 1am
h = t;
t = round(t);


fprintf('Select a number corresponding to the variable you''d like to calculate\n')
fprintf('A[0], B[1], C[2], k[3], E[4], delta[5], beta[6], phi(s)[7], IB[8], IBC[9], IDC[10], IRC[11], or IC[12]\n')
wh = 'what would you like to calculate: ?\n'; 
ww = input(wh);

if h >= 12
            H = -15*(h-12);
        else
            H = 15*(12-h);
end  
beta = asind(cosd(L)*cosd(de)*cosd(H)+sind(L)*sind(de));
if beta<=0
    beta = 0;
end
phis1 = asind((cosd(de)*sind(H))./(cosd(beta)));
phis2 = 180 - asind((cosd(de)*sind(H))./(cosd(beta)));
xx = cosd(H);
xx2 = tand(de)/tand(L);
if xx>=xx2
    phis=phis1;
else
    phis=phis2;
end
%m = 1/sind(beta);
m = (((708*sind(beta))^2)+1417)^.5-(708*sind(beta)); %<--corrected
tetha = acosd(cosd(beta)*cosd(phis-phic(1,1))*sind(sigma) + sind(beta)*cosd(sigma));
IB = A*exp(-k*m);
IBC = IB*cosd(tetha);
IDC = IB*C*(.5+cosd(sigma)/2);
IRC = rho(1,1)*IB*(sind(beta)+C)*(.5-cosd(sigma)/2); %<--corrected
IC = IBC + IDC + IRC;
if ww == 0
    fprintf('A: %f [W/m^2]\n', A)
elseif ww == 1
    fprintf('B: %f [N/A]\n', B)
elseif ww == 2
    fprintf('C: %f [N/A]\n', C)
elseif ww == 3 
    fprintf('k: %f [N/A]\n', k)    
elseif ww == 4                          
    fprintf('E: %f [minutes]\n', E)  
elseif ww == 5
    fprintf('delta: %f [degree]\n', de)
elseif ww == 6  
    fprintf('beta: %f [degree]\n', beta)
elseif ww == 7
    fprintf('phis: %f [degree]\n', phis)
elseif ww == 8
    fprintf('IB: %f [W/m^2]\n', IB)
elseif ww == 9
    fprintf('IBC: %f [W/m^2]]\n', IBC)
elseif ww == 10
    fprintf('IDC: %f [W/m^2]\n', IDC)
elseif ww == 11
    fprintf('IRC: %f [W/m^2]\n', IRC)
elseif ww == 12
    fprintf('IC: %f [W/m^2]\n', IC)
else 
fprintf('make a selection between 0 and 12, thanks\n');
end


