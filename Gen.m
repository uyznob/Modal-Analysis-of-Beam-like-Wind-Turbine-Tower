% This function generates elements as 
% input number of segments
function ge = Gen(n,data)
    %Size of data
    [r,c] = size(data);    
    %Number of divisions for first element
    m = round(0.2*n);
    if m == 0
        m=1;
    end
    %Define ge size
    ge = zeros((r-1)*n+m,c);
    %========================
    %Data for first segment
    %========================
    for i=1:m               
        %Thickness stays the same
        ge(i,3) = data(1,3);       
        %Height        
        ge(i,4) = data(1,4)/m;           
        %Start diameter
        ds = data(1,1);
        de = data(1,2);
        ge(i,1) = de + (ds-de)*(m-i+1)/m;         
    end
    %Elevation        
        ge(1,5) = 0;
        if m>=2
            for i=2:m
                ge(i,5) = ge(i-1,5) + ge(i-1,4);
            end      
        end
    %========================
    %Data for other elements
    %========================    
    for i=2:r               
        %Thickness stays the same
        for i3 = 1:n
            ge((i-1)*n + i3 - (n-m),3) = data(i,3);
        end
        %Height
        for i4 = 1:n
            ge((i-1)*n + i4 - (n-m),4) = data(i,4)/n;
        end
        %Elevation
        for i5 = 1:n
            ge((i-1)*n + i5 - (n-m),5) = ge((i-1)*n + i5 - (n-m)-1,5) ... 
                                        + ge((i-1)*n + i5 - (n-m)-1,4);
        end
        %Start diameter
        ds = data(i,1);
        de = data(i,2);
        for i1 = 1:n
            ge((i-1)*n + i1 - (n-m),1) = de + (ds-de)*(n-i1+1)/n;            
        end            
    end
    %End diameter
    for i=1:(r-1)*n+m-1
        ge(i,2) = ge(i+1,1);
    end
    ge((r-1)*n+m,2) = data(r,2);    
end
