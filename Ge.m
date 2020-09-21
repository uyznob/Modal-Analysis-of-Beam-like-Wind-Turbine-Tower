% This function generates elements as 
% input number of segments
function ge = Ge(n,data)
    %Size of data
    [r,c] = size(data);
    %First line data doesn't change
    %Because element is small (550 mm)
    ge = zeros((r-1)*n+1,c);
    ge(1,:) = data(1,:);
    for i=2:r
        %Thickness stays the same
        for i3 = 1:n
            ge((i-1)*n + i3 - (n-1),3) = data(i,3);
        end
        %Height
        for i4 = 1:n
            ge((i-1)*n + i4 - (n-1),4) = data(i,4)/n;
        end
        %Elevation
        for i5 = 1:n
            ge((i-1)*n + i5 - (n-1),5) = ge((i-1)*n + i5 - (n-1)-1,5) ... 
                                        + ge((i-1)*n + i5 - (n-1)-1,4);
        end
        %Start diameter
        ds = data(i,1);
        de = data(i,2);
        for i1 = 1:n
            ge((i-1)*n + i1 - (n-1),1) = de + (ds-de)*(n-i1+1)/n;            
        end            
    end
    %End diameter
    for i=2:(r-1)*n
        ge(i,2) = ge(i+1,1);
    end
    ge((r-1)*n+1,2) = data(r,2);    
end
