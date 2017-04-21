 function[sulfonate_owner]=sulfonate_ownership_function(x_H3O, y_H3O, z_H3O,  x_sulf,y_sulf,z_sulf,...
     Radius_si, x_length_si,scale,r_trap)

%% Initialize
% Radius_si=1;
% x_length_si=10;

%Values from nanometers to ticks
Radius=round(Radius_si*scale);
radius_min_2=((Radius_si-r_trap)*scale-1)^2;%Protons outside this radius are in the surface
% x_length=round(x_length_si*scale);
circ=2*pi*Radius;

rad_2=y_H3O(:).^2+z_H3O(:).^2;
owned=length(find(rad_2>radius_min_2));

sulfonate_owner=zeros(2,owned);
owner_index=1;

for i=1:length(x_H3O(:))
	if rad_2(i)>=radius_min_2
		d_low_2=circ^2; %Set an unrealistically high starting minimum distance
		index_close=1; %Index of the closest sulfonate
        for j=1:length(x_sulf)
			d_i_2=(x_H3O(i)-x_sulf(j))^2+(y_H3O(i)-y_sulf(j))^2+(z_H3O(i)-z_sulf(1,j))^2;
            already=find(sulfonate_owner(2,:)==j, 1);
            if d_i_2<d_low_2 && isempty(already)
                d_low_2=d_i_2;
                index_close=j;
            end
        end
        sulfonate_owner(1,owner_index)=i;
        sulfonate_owner(2,owner_index)=index_close;
        owner_index=owner_index+1;
	end
end

return;