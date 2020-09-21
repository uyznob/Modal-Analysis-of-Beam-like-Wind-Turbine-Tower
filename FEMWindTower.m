clear all;
clc;
%===============================
%All units in N, m, kg
%Note the first 0.55 segment
%is not devided because. 0.6m 
%approaches the stable natural
%frequency
%===============================

filename = 'Finite Element Data.xlsx';
sheet = 1;
xlRange = 'B3:F30';
data = xlsread(filename,sheet,xlRange);
%Divide the segments into n smaller segments
datan = Gen(2,data);
% datan = round(datan);
[r,c] = size(datan);
%Define global stiffness and mass matrix
Kg = zeros((r+1)*2);
Mg = zeros((r+1)*2);
%Define elastic modulus N/m2
E = 2.1e+011; 
%Establish global stiffness and mass matrix
for i=1:r
    %In meter
    ds = datan(i,1)/1000;
    de = datan(i,2)/1000;
    t  = datan(i,3)/1000;
    h  = datan(i,4)/1000;
    
    %Avarage outer diameter m
    dao = (ds+de)/2+t;
    %Average inner diameter m
    dai = (dao - t*2);    
    %Moment of inertia m^4
    I = pi*((dao)^4-(dai)^4)/64;
    %Mass with density 7850 kg/m^3
    m = pi/4*((dao)^2-(dai)^2)*h*7850;
    %Uniform mass
    m = m/h;
    
    %Stiffness of element: Chopra P.741        
    ke = E*I/h^3.*[12 6*h -12 6*h; 6*h 4*h^2 -6*h 2*h^2; 
                  -12 -6*h 12 -6*h; 6*h 2*h^2 -6*h 4*h^2];              
    %Mass of element: Chopra P.742   kg/m
    me = m*h/420.*[156 22*h 54 -13*h; 22*h 4*h^2 13*h -3*h^2;
                54 13*h 156 -22*h; -13*h -3*h^2 -22*h 4*h^2];  
    %The following me is also working but only for dynamics
    %me = m*h.*[0.5 0 0 0;0 0 0 0;0 0 0.5 0;0 0 0 0];         
    
    %Global stiffness matrix
    ri = 2*i-1; ci=2*i+2;
    Kg(ri:ci,ri:ci) = ...
        Kg(ri:ci,ri:ci)+ke;    
    %Global mass matrix
    Mg(ri:ci,ri:ci) = ...
        Mg(ri:ci,ri:ci)+me;  
end
%Add lump masss including rotor and nacelle
Mg((r+1)*2-1,(r+1)*2-1) = Mg((r+1)*2-1,(r+1)*2-1) + 107.8e3;
Mg((r+1)*2-3,(r+1)*2-3) = Mg((r+1)*2-3,(r+1)*2-3) + 107.8e3;
%Add support condition
Kgn = Kg(3:(r+1)*2,3:(r+1)*2);
Mgn = Mg(3:(r+1)*2,3:(r+1)*2);
%Solve eigen value problem
[Psi, w2] = eig(Kgn,Mgn);
f = sqrt(diag(w2))/(2*pi);

%===============================
%Plot the result
%Note: Not plot the (0,0) point
%===============================
%Find elevation of element
z = zeros(r,1);
z(:,1) = datan(:,5)/1000;
z(r,1) = sum(datan(:,4))/1000;
%Extract horizontal displacement ux
ux = zeros(r,5);
for i=1:r
    ux(i,1)= Psi(2*(i-1)+1,1);
    ux(i,2)= Psi(2*(i-1)+1,2);
    ux(i,3)= Psi(2*(i-1)+1,3);
    ux(i,4)= Psi(2*(i-1)+1,4);
    ux(i,5)= Psi(2*(i-1)+1,5);
end
%Plot first 5 mode shapes
figure
subplot(1,5,1)
plot(ux(:,1),z)
hold on
plot(0,z,'m-.')
subplot(1,5,2)
plot(ux(:,2),z)
hold on
plot(0,z,'m-.')
subplot(1,5,3)
plot(ux(:,3),z)
hold on
plot(0,z,'m-.')
subplot(1,5,4)
plot(ux(:,4),z)
hold on
plot(0,z,'m-.')
subplot(1,5,5)
plot(ux(:,5),z)
hold on
plot(0,z,'m-.')
%Show first 5 natural frequency
disp('First 5 natural frequency:');
disp(f(1:5,1));