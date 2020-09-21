%The location of flange
function fla = fla(n)      
    %Number of divisions for first element
    m = round(0.2*n);
    if m == 0
        m=1;
    end
    %Define fla size
    fla = zeros(7,4);
    %Element ID contains 1st flange
    fla(1,1) = m;       %End of this is flange
    fla(2,1) = m+1;     %Beginning of this is flange
    %Other element ID contains remaining flanges
    fla(3,1) = m+7*n;
    fla(4,1) = m+7*n+1;
    fla(5,1) = m+17*n;
    fla(6,1) = m+17*n+1;
    fla(7,1) = m+27*n;
    %Flange information
    fla(:,2:4) =  [ 90 230 38;
                   130 250 40;
                   130 250 26;
                   130 250 26;
                   110 220 21;
                   110 220 21;
                   250 123 16];    
end
