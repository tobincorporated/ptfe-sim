function [proceed]=sulfonate_transfer(x,y,z,x_sulf,y_sulf,z_sulf, Radius, scale,sulfonate_owner_i)

%Purpose of this function is to determine if a desired surface 
%positive charge can move as desired. Must work for Grotthuss and 
%translation


%Algorithm:
%1. What sulfonate does it belong to?
%2. Find where the hydronium wants to move
%3. Is this closer to a different sulfonate?
%4. Is this outside the "surface", in the bulk?
%5. If 3 && 4 ==no, whatever, it can do it. 
%6. If it wants to move to a different sulfonate, rand<e^En/kT
%7. Then it can move.
%8. Return value:  proceed==0,1

%%

proceed=1;

% if z^2+y^2 <= (Radius-scale*1)^2
% 	if rand < 1
% 		proceed=0;
% 	end
% end

kbe=10^-23;
kb=1.38*kbe; %Boltzman's constant, J/K
T=300; %Temperature, K
ele=10^-19;
el=1.602*ele; %Electron charge, C
ep_e=10^-12;
ep_0 = 8.854*ep_e; %Electric constant epsilon_0, C^2/N*m^2 
Re=10^-9;% scale of nanometer
Rf=.244 *Re; %Radius of sulfonate
Ri=.143 *Re; %Radius of H3O
li=.255*Re; %Distance of energy barrier
kappa=6;

En0=-el^2/(4*pi*ep_0*kappa) *1/(Rf+Ri);
En1=-el^2/(4*pi*ep_0*kappa) *1/(Rf+Ri+li);

exp0=exp(-En0/(kb*T));
exp1=exp(-En1/(kb*T));

% P1=exp1/(exp0);
P1=1;

d_0_2=(x-x_sulf(sulfonate_owner_i))^2+(y-y_sulf(sulfonate_owner_i))^2+(z-z_sulf(sulfonate_owner_i))^2;

d_low=d_0_2;

for i=1:length(x_sulf)
	d_1_2=(x-x_sulf(i))^2+(y-y_sulf(i))^2+(z-z_sulf(i))^2;
    if d_1_2 < d_low
		if rand < P1
			d_low=d_1_2;
            proceed=1;
		else
			proceed=0;
		end
    end		
end

return;