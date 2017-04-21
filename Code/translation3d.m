function  [x_H2Omas, y_H2Omas, z_H2Omas, x_H3Omas, y_H3Omas, z_H3Omas, mid_current,sulf_current]=...
    translation3d(x_H2Omas, y_H2Omas, z_H2Omas, x_H3Omas, y_H3Omas, z_H3Omas, x_sulf, y_sulf, z_sulf,...
    N_x, Radius_nm ,sulfonate_owner,H2O_grid, E_grid,r_trap,exp_barrier,mode)

scale=10/.5774;
mid_current=0; %Initialize current values
sulf_current=0;

Radius=Radius_nm*scale; %Radius in ticks
radmax=floor(Radius);

for i=1:length(x_H2Omas(1,:))

%% Set H2O direction
            x_0=x_H2Omas(i);
            y_0=y_H2Omas(i);
            z_0=z_H2Omas(i);

            proceed=1;
            T_dir = ceil(rand*6);
            x_plus=0;
            y_plus=0;
            z_plus=0;
            
            
            if x_0<1
                   T_dir=1;
            end

            if x_0>N_x
                   T_dir=2;
            end
            
            switch T_dir
                case 1
                  x_plus=1;
                case 2
                  x_plus=-1;
                case 3
                  y_plus=1;
                case 4
                  y_plus=-1;  
                case 5
                  z_plus=1;
                case 6
                  z_plus=-1;                              
            end

            
             x=x_0+x_plus;
             y=y_0+y_plus;
             z=z_0+z_plus;
             
%% Avoid boundaries


            if (y)^2+(z)^2>=(radmax-1)^2
                proceed=0;
            end
             
%% Avoid collisions
        
        if proceed==1
            %O_spot adjusts the position to a grid location, plus a short
            %buffer to keep molecules away from each other
            O_spot_x=x+3;
            O_spot_y=y+radmax+3;
            O_spot_z=z+radmax+3;
            
            O_spot=H2O_grid(O_spot_x,O_spot_y,O_spot_z);

            %Because O_spot is 2 steps away from the oxygen, it is guaranteed
            %not to overlap the original surroundings of the molecule
            if O_spot>0
                proceed=0;            
            end
        end
      
        
%% set values and reloop
        if proceed==1
            x_H2Omas(i)= x_H2Omas(i)+ x_plus ;
            y_H2Omas(i)= y_H2Omas(i)+ y_plus ;
            z_H2Omas(i)= z_H2Omas(i)+ z_plus ;
        end
end

%% H3O ions

for i=1:length(x_H3Omas(:))

            x_0=x_H3Omas(1,i);
            y_0=y_H3Omas(1,i);
            z_0=z_H3Omas(1,i);

            proceed=1;
%% Choose direction to translate
            
            
            x_plus=0;
            y_plus=0;
            z_plus=0;
            
            T_dir=ceil(6*rand);

            
            switch T_dir
                case 1
                    x_plus=1;
                case 2
                    x_plus=-1;
                case 3
                    y_plus=1;
                case 4
                    y_plus=-1;
                case 5 
                    z_plus=1;
                case 6 
                    z_plus=-1;
            end
%% Avoid boundaries

        if (y_0)^2+(z_0)^2 >= (radmax-2)^2 
            if abs(y_0)>abs(z_0)
                if y_0>0
                    y_plus=-1;
                    z_plus=0;
                    x_plus=0;
                else
                    y_plus=1;
                    z_plus=0;
                    x_plus=0;
                end
            else if z_0>0
                    z_plus=-1;
                    y_plus=0;
                    x_plus=0;
                 else
                    z_plus=1;
                    y_plus=0;
                    x_plus=0;
                end
            end
            proceed=1;
        end
    
        
            if x_0<2 && x_plus ==-1 
               x_plus=N_x-2;
            end

            if x_0>N_x-1 && x_plus==1
              x_plus=2-N_x;
            end
            
            x=x_0+x_plus;
            y=y_0+y_plus;
            z=z_0+z_plus;
            

%% Monte Carlo for energies and distribution. 
                
         

            E1=E_grid(x+3,y+radmax+3,z+radmax+3);
            E0=E_grid(x_0+3,y_0+radmax+3,z_0+radmax+3);
            

                if x_plus>=1
                    E0=E_grid(3,y+radmax+3,z+radmax+3);
                    E1=E_grid(2,y+radmax+3,z+radmax+3);
                end

                if x_plus<=-1
                    E0=E_grid(2,y+radmax+3,z+radmax+3);
                    E1=E_grid(3,y+radmax+3,z+radmax+3);
                end
 
            
            prob=exp(E0-E1);
            rand_flux=rand;
            
            if rand_flux > prob
                 proceed=0;
            end
%             
%             if x_plus ~=0
%                 pause;
%             end
%             



%% Avoid collisions
    
        %O_spot adjusts the position to a grid location, plus a short
        %buffer to keep molecules away from each other
        x_proxy=0;
        if abs(x_plus)==1
            x_proxy=x_plus/abs(x_plus);
        elseif x_plus>1
            x_proxy=-1;
        elseif x_plus<-1
            x_proxy=1;
        end
        O_spot_x=x+3+x_proxy;
        O_spot_y=y+radmax+3+y_plus;
        O_spot_z=z+radmax+3+z_plus;
if O_spot_x<1 || O_spot_y<1 ||O_spot_z<1 
    pause;
end
        O_spot=H2O_grid(O_spot_x,O_spot_y,O_spot_z);
        
        %Because O_spot is 2 steps away from the oxygen, it is guaranteed
        %not to overlap the original surroundings of the molecule
        if O_spot>0
            proceed=0;     
        end
        
        
%% Sulfonate hopping

%Energy barrier to leave owner sulfonate

    if proceed==1 && (y^2+z^2>(Radius-r_trap*scale-2)^2  )

        own_index=find(sulfonate_owner(1,:)==i, 1);
        if ~isempty(own_index)             
            sulf_i=sulfonate_owner(2,own_index);
            x_sulf_i=x_sulf(sulf_i);
            y_sulf_i=y_sulf(sulf_i);
            z_sulf_i=z_sulf(sulf_i);
            
            r0_2=(x_0-x_sulf_i)^2+(y_0-y_sulf_i)^2+(z_0-z_sulf_i)^2;
            r_2=(x-x_sulf_i)^2+(y-y_sulf_i)^2+(z-z_sulf_i)^2;
            if r_2>(r_trap*scale)^2 && r0_2<(r_trap*scale)^2 && rand>exp(exp_barrier)%% Energy barrier
                proceed=0;
            end
            
            if proceed==1
               if x-x_0~=0
                   x_disp=x-x_0;
                   sulf_current=sulf_current+x_disp/abs(x_disp); 
               end
            end
            
        end

    end
        
        
        
        
%% set values and reloop
        if proceed==1
            if abs(x_plus)==1
                mid_current=mid_current+x_plus;
            elseif abs(x_plus)>1
                mid_current=mid_current-x_plus/abs(x_plus);
            
            end
            
            x_H3Omas(i)= x_H3Omas(i)+ x_plus ;
            y_H3Omas(i)= y_H3Omas(i)+ y_plus ;
            z_H3Omas(i)= z_H3Omas(i)+ z_plus ;
        end

    
    
                 
                 
end


mid_current=mid_current/N_x;%%Average the current
sulf_current=sulf_current/N_x;

return;




%% 
E2grid=zeros(length(E_grid(1,:,1)));
sumgrid=sum(E_grid);
for i=1:74
    for j=1:74
        
        E2grid(i,j)=abs(sumgrid(1,i,j));
    end
end
surf(E2grid)