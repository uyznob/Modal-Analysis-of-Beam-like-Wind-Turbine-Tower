clear all;
clc;
%===============================
%All units in N, m, kg
%Note the first 0.55 segment
%===============================

filename = 'Finite Element Data.xlsx';
sheet = 1;
xlRange = 'B3:F30';
data = xlsread(filename,sheet,xlRange);
%Divide the segments into n smaller segments
datan = Gen(80,data);
datan = round(datan);
% datan = round(datan);
[r,c] = size(datan);
%Define global stiffness and mass matrix
Kg = zeros((r+1)*3);
Mg = zeros((r+1)*3);
%Define elastic modulus N/m2
E = 2.1e+011; 
%Establish global stiffness and mass matrix
ke = zeros(6,6);
me = zeros(6,6);
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
    %Area of section
    A = pi*((dao)^2-(dai)^2)/4;
    %Moment of inertia m^4
    I = pi*((dao)^4-(dai)^4)/64;
    %Volum of segment
    Vout= pi*h/3*((0.5*(ds+t))^2+(0.5*(de+t))^2+0.25*(ds+t)*(de+t));
    Vin = pi*h/3*((0.5*(ds-t))^2+(0.5*(de-t))^2+0.25*(ds-t)*(de-t));
    V   = Vout-Vin;
    %Mass with density 7698 kg/m^3
    m = V*7698;
    %Uniform mass
    m = m/h;
    
    %Stiffness of element: Chopra P.741        
    k1  = E*A/h.*[1 -1; -1 1]; 
    k23 = E*I/h^3.*[12 6*h -12 6*h; 6*h 4*h^2 -6*h 2*h^2; 
                  -12 -6*h 12 -6*h; 6*h 2*h^2 -6*h 4*h^2];
    %Assemble to ke
    ke(1,1) = k1(1,1);ke(1,4) = k1(1,2);
    ke(4,1) = k1(2,1);ke(4,4) = k1(2,2);
    
    ke(2,5) = k23(1,3);ke(2,6) = k23(1,4);
    ke(3,5) = k23(2,3);ke(3,6) = k23(2,4);
    
    ke(5,2) = k23(3,1);ke(6,2) = k23(4,1);
    ke(5,3) = k23(3,2);ke(6,3) = k23(4,2);
    for ir= 1:2
        for ic= 1:2
            ke(ir+1,ic+1) = k23(ir,ic);
        end
    end
    for ir= 3:4
        for ic= 3:4
            ke(ir+2,ic+2) = k23(ir,ic);
        end
    end                   
    %Mass of element: Chopra P.742 (kg/m)
    m1  = m*h/420.*[140 70;70 140];
    m23 = m*h/420.*[156 22*h 54 -13*h; 22*h 4*h^2 13*h -3*h^2;
                54 13*h 156 -22*h; -13*h -3*h^2 -22*h 4*h^2];
    %Assemble to me
    me(1,1) = m1(1,1);me(1,4) = m1(1,2);
    me(4,1) = m1(2,1);me(4,4) = m1(2,2);
    
    me(2,5) = m23(1,3);me(2,6) = m23(1,4);
    me(3,5) = m23(2,3);me(3,6) = m23(2,4);
    
    me(5,2) = m23(3,1);me(6,2) = m23(4,1);
    me(5,3) = m23(3,2);me(6,3) = m23(4,2);
    for ir= 1:2
        for ic= 1:2
            me(ir+1,ic+1) = m23(ir,ic);
        end
    end
    for ir= 3:4
        for ic= 3:4
            me(ir+2,ic+2) = m23(ir,ic);
        end
    end     
    %Global stiffness matrix
    ri = 3*i-2; ci=3*i+3;
    Kg(ri:ci,ri:ci) = ...
        Kg(ri:ci,ri:ci)+ke;    
    %Global mass matrix
    Mg(ri:ci,ri:ci) = ...
        Mg(ri:ci,ri:ci)+me;  
end
%Add lump masss including rotor and nacelle
Mg((r+1)*3-1,(r+1)*3-1) = Mg((r+1)*3-1,(r+1)*3-1) + 107.8e3;
Mg((r+1)*3-2,(r+1)*3-2) = Mg((r+1)*3-2,(r+1)*3-2) + 107.8e3;
Mg((r+1)*3-4,(r+1)*3-4) = Mg((r+1)*3-4,(r+1)*3-4) + 107.8e3;
Mg((r+1)*3-5,(r+1)*3-5) = Mg((r+1)*3-5,(r+1)*3-5) + 107.8e3;
%Add support condition
Kgn = Kg(4:(r+1)*3,4:(r+1)*3);
Mgn = Mg(4:(r+1)*3,4:(r+1)*3);
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
ux = zeros(r,3);
for i=1:r
    ux(i,1)= Psi(3*(i-1)+2,1);
    ux(i,2)= Psi(3*(i-1)+2,2);
    ux(i,3)= Psi(3*(i-1)+2,3);
    ux(i,4)= Psi(3*(i-1)+2,4);
    ux(i,5)= Psi(3*(i-1)+2,5);
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