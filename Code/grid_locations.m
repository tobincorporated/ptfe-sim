function [H2O_grid]=grid_locations(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, x_max,y_max,z_max)

H2O_grid=zeros(x_max+6, 2*y_max+6,2*z_max+6);

%Grid locations as a decimal number.  The units correspond to the number of
%the molecule in the array, the decimal corresponds to the atom type
%    H2O
% .1 = H2O space
% .2 = H2O oxygen
% .3 = H2O hydrogen
%    H3O
% .4 = H3O space
% .5 = H3O oyxgen
% .6 = H3O hydrogen
% eg. 12.3  means molecule 12, atom type 3 (h3o oxygen atom)


for i=1:length(x_H2O(:))
    x_O=x_H2O(i)+3;
    y_O=y_H2O(i)+y_max+3;
    z_O=z_H2O(i)+z_max+3;
    
    for ii=-1:1
        for jj=-1:1
            for kk=-1:1
                H2O_grid(x_O+ii, y_O+jj,z_O+kk)=i+.1;
            end
        end
    end
    
    H2O_grid(x_O,y_O,z_O)= i+.2;
   
    
    

   
end

for i=1:length(x_H3O)

   
    
    x_O=x_H3O(i)+3;
    y_O=y_H3O(i)+y_max+3;
    z_O=z_H3O(i)+z_max+3;
    
    
    for ii=-1:1
        for jj=-1:1
            for kk=-1:1
                H2O_grid(x_O+ii, y_O+jj,z_O+kk)=i+.4;
            end
        end
    end
    
    
    H2O_grid(x_O,y_O,z_O)= i+.5;
    
 
   
end


return;