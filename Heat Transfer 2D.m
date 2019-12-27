clc
close all 
clear all

k = 52.35; %W/m-K
C = 523; %J/kg-K
dx = 0.0254;%m
rho = 7854; %kg/m^3
dtmin = rho*C*dx^2/4/k %s Calculated from Stability of node 3

dt = 12.65561 %s User input
Fo = k*dt/rho/C/dx^2
ptotal = round(3600/dt);
ttotal = dt*ptotal/60; %min

h1 = 31; %W/m^2-K
Bi1 = h1*dx/k 

h2 = 17.2; %W/m^2-K
Bi2 = h2*dx/k

Ti = 500; %oK
Tinf1 = 400; %oK
Tinf2 = 320; %oK

A = zeros(20)
Tp = zeros(20,1)

Tp(1:18) = Ti;
Tp(19) = Tinf1;
Tp(20) = Tinf2;

%% Tinf
A(19,19) = 1;
A(20,20) = 1;

%% Node 1
A(1,1) = 1-4*Bi1*Fo-2*Fo;
A(1,3) = 2*Fo;
A(1,19) = 4*Bi1*Fo;
A(1,20) = 0;

%% Node 2
A(2,2) = 1-2*Fo-Bi1*Fo-Bi2*Fo;
A(2,3) = 2*Fo;
A(2,19) = Bi1*Fo;
A(2,20) = Bi2*Fo;

%% Node 3
A(3,3) = 1-4*Fo;
A(3,1) = Fo;
A(3,2) = Fo;
A(3,4) = Fo;
A(3,5) = Fo;

%% Node 4
A(4,4) = 1-2*Fo-Bi1*Fo-Bi2*Fo;
A(4,3) = 2*Fo;
A(4,19) = Bi1*Fo;
A(4,20) = Bi2*Fo;

%% Node 5
A(5,5) = 1-2*Fo-2*Bi2*Fo;
A(5,3) = Fo;
A(5,6) = Fo;
A(5,20) = 2*Bi2*Fo;

%% Node 6-17
A(6,6) = 1-2*Fo-2*Bi2*Fo;
A(6,5) = Fo;
A(6,7) = Fo;
A(6,20) = 2*Bi2*Fo;

%% Node 18
A(18,18) = 1-2*Fo-2*Bi2*Fo;
A(18,17) = 2*Fo;
A(18,20) = 2*Bi2*Fo;

%% Node 6-17
for i = 6:17
    A(i,i-1)= Fo;
    A(i,i)= 1-2*Fo-2*Bi2*Fo;
    A(i,i+1)= Fo;
    A(i,20) = 2*Bi2*Fo;
end

 %% 2D Matrix
B = [Tp(19) Tp(1) Tp(1) Tp(19); Tp(19) Tp(1) Tp(1) Tp(19); Tp(2) Tp(3) Tp(3) Tp(4); Tp(2) Tp(3) Tp(3) Tp(4); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(6) Tp(6) Tp(20);
    Tp(20) Tp(6) Tp(6) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(9) Tp(9) Tp(20);
    Tp(20) Tp(9) Tp(9) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(12) Tp(12) Tp(20);
    Tp(20) Tp(12) Tp(12) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(15) Tp(15) Tp(20);
    Tp(20) Tp(15) Tp(15) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(18) Tp(18) Tp(20)]

%% Plot and Capture GIF
X = [0, 0.5, 1.5, 2];
Y = 15:-0.5:0;
v = 300:20:500;
contourf(X,Y, [Tp(19) Tp(1) Tp(1) Tp(19); Tp(19) Tp(1) Tp(1) Tp(19); Tp(2) Tp(3) Tp(3) Tp(4); Tp(2) Tp(3) Tp(3) Tp(4); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(6) Tp(6) Tp(20);
    Tp(20) Tp(6) Tp(6) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(9) Tp(9) Tp(20);
    Tp(20) Tp(9) Tp(9) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(12) Tp(12) Tp(20);
    Tp(20) Tp(12) Tp(12) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(15) Tp(15) Tp(20);
    Tp(20) Tp(15) Tp(15) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(18) Tp(18) Tp(20)], v);

title(['Temp(oK), n=0, t=0s'])
ylabel('Height (in)'), xlabel('Width(in)')
colorbar; 
frame = getframe(1)
im{1} = frame2im(frame)
colormap(jet(500))

for j = 1:ptotal;
    hold all    
    Tp1 = A*Tp;
    title(['Temp(^oK), n=' num2str(j) ', t=' num2str(dt*j) 's'])
    h = contourf(X,Y, [Tp(19) Tp(1) Tp(1) Tp(19); Tp(19) Tp(1) Tp(1) Tp(19); Tp(2) Tp(3) Tp(3) Tp(4); Tp(2) Tp(3) Tp(3) Tp(4); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(5) Tp(5) Tp(20); Tp(20) Tp(6) Tp(6) Tp(20);
    Tp(20) Tp(6) Tp(6) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(7) Tp(7) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(8) Tp(8) Tp(20); Tp(20) Tp(9) Tp(9) Tp(20);
    Tp(20) Tp(9) Tp(9) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(10) Tp(10) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(11) Tp(11) Tp(20); Tp(20) Tp(12) Tp(12) Tp(20);
    Tp(20) Tp(12) Tp(12) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(13) Tp(13) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(14) Tp(14) Tp(20); Tp(20) Tp(15) Tp(15) Tp(20);
    Tp(20) Tp(15) Tp(15) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(16) Tp(16) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(17) Tp(17) Tp(20); Tp(20) Tp(18) Tp(18) Tp(20)],v);   
    drawnow
    frame = getframe(1);
    im{j+1} = frame2im(frame);   
    Tp = Tp1;  
end

filename = 'TransientPlaneWall.gif'; % Specify the output file name
for k = 1:j+1
    [A,map] = rgb2ind(im{k},256);
    if k == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.05);
    end
end
