function  [x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, mid_current, sulf_current]=...
    grotthuss3d(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, x_sulf, y_sulf, z_sulf,... 
    N_x, Radius_si, P_G, sulfonate_owner,H2O_grid, E_grid, r_trap,exp_barrier)


scale=10/.5774;
Radius=Radius_si*scale;
radmax=floor(Radius);
mid_current=0;
sulf_current=0;
hoplength_nm=.3;%Max distance from ion to molecule for Grotthuss (nm)
hoplength=hoplength_nm*scale;%In tick marks

for i=1:length(x_H3O(:))
    if rand<P_G
        proceed=1;
        H2O_indeces=zeros(30);
        H2O_index=0;
       
            %look in nearby space for H2O         
            x=x_H3O(i)+3;
            y=y_H3O(i)+3+radmax;
            z=z_H3O(i)+3+radmax;
            
            %For loops over to check the corners of the H3O molecule's box
            %space for overlap with an H2O molecule
            for ii=-hoplength:hoplength:hoplength
                for jj=-hoplength:hoplength:hoplength
                    for kk=-hoplength:hoplength:hoplength
                        x1=round(x+ii);
                        y1=round(y+jj);
                        z1=round(z+kk);
                        
                        if x1>N_x+6
                            x1=N_x+6;
                        end
                        if x1<1
                            x1=1;
                        end
                        if y1>2*radmax+6
                            y1=2*radmax+6;
                        end
                        if y1<1
                            y1=1;
                        end
                        if z1>2*radmax+6
                            z1=2*radmax+6;
                        end
                        if z1<1
                            z1=1;
                        end
                        
                        H2O_num=H2O_grid(x1,y1,z1);
                        H2O_num_floor=floor(H2O_num);%# of the molecule
                        H2O_num_type=mod(H2O_num,1);% type of molecule
                        
                        %1)If the space is occupied
                        %2)If it is H2O
                        if H2O_num>0 && H2O_num_type<.39 
                            
                            %grid location of h2o
                            xH1=x_H2O(H2O_num_floor)+3;
                            yH1=y_H2O(H2O_num_floor)+3+radmax;
                            zH1=z_H2O(H2O_num_floor)+3+radmax;
                            
                            d_1_2=(x-xH1)^2+(y-yH1)^2+(z-zH1)^2;
                            
                            if d_1_2<hoplength^2
                                H2O_index=H2O_index+1;
                                H2O_indeces(H2O_index)=H2O_num_floor;
                                %catalogs the close enough molecules
                            end
                            
                        end

                    end
                end
            end
        
    if H2O_index==0
        proceed=0;
    else
        H2O_index_spot=ceil(rand*H2O_index);%Choose randomly which h2o
        H2Onumber=H2O_indeces(H2O_index_spot);%picks out the h2o number to use
        
    end
    
    if proceed==1

        x0=x_H3O(i);
        y0=y_H3O(i);
        z0=z_H3O(i);

        xH=x_H2O(H2Onumber);
        yH=y_H2O(H2Onumber);
        zH=z_H2O(H2Onumber);
    end

%% Sulfonate hopping

%Energy barrier to leave owner sulfonate

    if proceed==1 && (y0^2+z0^2>(Radius-r_trap*scale-2)^2  )

        
        own_index=find(sulfonate_owner(1,:)==i, 1);
        if ~isempty(own_index)             
            sulf_i=sulfonate_owner(2,own_index);
            
            x_sulf_i=x_sulf(sulf_i);
            y_sulf_i=y_sulf(sulf_i);
            z_sulf_i=z_sulf(sulf_i);
            
            r0_2=(x0-x_sulf_i)^2+(y0-y_sulf_i)^2+(z0-z_sulf_i)^2;
            r_2=(xH-x_sulf_i)^2+(yH-y_sulf_i)^2+(zH-z_sulf_i)^2;
            if r_2>(r_trap*scale)^2 && r0_2<(r_trap*scale)^2 && rand>exp(exp_barrier)%% Energy barrier
                proceed=0;
            end
%% Around here
            if proceed==1 && r_2>(r_trap*scale)^2 && r0_2<(r_trap*scale)^2
               if xH-x0~=0
                   x_disp=xH-x;
                   sulf_current=sulf_current+x_disp/abs(x_disp); 
               end
            end
        
            
        end

    end


%% Energies


         

            if proceed==1
                if z0> radmax
                    pause(.1);
                end
                E1=E_grid(xH+3,yH+radmax+3,zH+radmax+3);
                E0=E_grid(x0+3,y0+radmax+3,z0+radmax+3);

    %             E0=-.15*(y0^2+z0^2).^2/(scale^4);
    %             E1=-.15*(y^2+z^2).^2/(scale^4);

                prob=exp(E0-E1);
    %             prob=1;
                rand_flux=rand;

                if rand_flux > prob
                     proceed=0;
                end
            end
%% If proceed==1, let the Grotthuss mechanism occur


            if proceed==1

                
                x_2=x_H2O(H2Onumber);
                x_3=x_H3O(i);
                
                hop_disp=(x_2-x_3);
                if hop_disp==0
                    current_plus=0;
                else current_plus=hop_disp;
                    
                end
                
                mid_current=mid_current+current_plus;

                x_temp=x_H3O(i);
                y_temp=y_H3O(i);
                z_temp=z_H3O(i);

                x_H3O(i)=x_H2O(H2Onumber);
                y_H3O(i)=y_H2O(H2Onumber);
                z_H3O(i)=z_H2O(H2Onumber);

                x_H2O(H2Onumber)=x_temp;
                y_H2O(H2Onumber)=y_temp;
                z_H2O(H2Onumber)=z_temp;
            end  
        
    end
end


mid_current=mid_current/N_x;
sulf_current=sulf_current/N_x;
return;