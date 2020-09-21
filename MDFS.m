clear all;
clc;
%===============================
%All units in N, m, kg
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
Kg = zeros(r);
Mg = zeros(r);
%Define elastic modulus N/m2
E = 2.1e+011; 
%Establish global stiffness and mass matrix
for i=2:r
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
    %Stiffness of element: Chopra P.741        
    ke = 12*E*I/h^3.*[1 -1;-1 1];    
    %Global stiffness matrix    
    Kg(i-1:i,i-1:i) = ...
        Kg(i-1:i,i-1:i)+ke;    
    %Global mass matrix
    Mg(i,i) = m;  
end
ds = datan(1,1)/1000;
de = datan(1,2)/1000;
t  = datan(1,3)/1000;
h  = datan(1,4)/1000;

%Avarage outer diameter m
dao = (ds+de)/2+t;
%Average inner diameter m
dai = (dao - t*2);    
%Moment of inertia m^4
I = pi*((dao)^4-(dai)^4)/64;
%Mass with density 7850 kg/m^3
k1 = 12*E*I/h^3;
m1 = pi/4*((dao)^2-(dai)^2)*h*7850;  
Kg(1,1) = Kg(1,1) + k1;
Mg(1,1) = m1;
%Add lump masss including rotor and nacelle
Mg(r,r) = Mg(r,r) + 107.8e3;

%Solve eigen value problem
[Psi, w2] = eig(Kg,Mg);
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
    ux(i,1)= Psi(i,1);
    ux(i,2)= Psi(i,2);
    ux(i,3)= Psi(i,3);
    ux(i,4)= Psi(i,4);
    ux(i,5)= Psi(i,5);
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