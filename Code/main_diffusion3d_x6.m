function [E_g , exp_barrier,  Radius , E_strength,  current_tot  ,...
    current_d ,  current_g , del_rho, current_ticks, current_density,...
    lambda, N, rho_ave, steady_rho_rad]= main_diffusion3d_x6(Radius, E_strength, mode)

%Initialization
global t;
tic;
rand('twister', sum(100*clock)); % Initialize random functions


%Parameters
%length in nanometers
x_length=3;%length of the channel
t_max=6002; % Ensures the program doesn't loop endlessly 
current_ticks=1000*(floor(t_max/1000)-1);% timesteps over which current is recorded

r_trap=.4; %Trap radius for sulfonates
E_g= 2;% Strength of Grotthuss Mechanism
P_G=exp(-E_g);
exp_barrier=-.1;% Sulfonate bond barrier in kbT
target_J=1/8/(pi*Radius^2);%Target current density per timestep, ions/nm^2/step

target_I=target_J*pi*Radius^2;

%Grid is scaled so that a molecular bond from oxygen to hydrogen
%is the three-dimensional diagonal of a 1x1x1 cube in the grid.  
%This length is approximately 1 Angstrom

%1 mark= sqrt(1/3) angstrom ~ .5774 angstroms
%10 angstrom = 1 nanometer
%10 marks = sqrt(1/3) nanometer
scale=10/.5774;  %scale = ticks/nm
timescale=1;%timescale in seconds. Determine with kinetic theory


%Maximum and minimum tick mark values, converted from nanometers to ticks

radmax=floor(Radius*scale);
radmaxmax=floor(3*scale);
x_max=floor(x_length*scale);
y_max=radmax;
z_max=radmax;

t=1; % t= timestep #left_current=zeros(1,t_max); %Create empty vector for the current through the left border

mid_currentd=zeros(1,t_max); %Create empty vector for the current through the middle, diffusive
mid_currentg=zeros(1,t_max); %As above, but for Grotthuss
mid_currentb=zeros(1,t_max); %As above, but for borders
sulf_currentd=zeros(1,t_max); %diffusive  sulfonate current
sulf_currentg=zeros(1,t_max); %grotthuss  sulfonate current

[x_sulf, y_sulf, z_sulf]=sulfonate_place(Radius, x_length, scale);


%% Number of molecules

%Average of 33.4 molecules per cubic nanometer, with Gaussian random
%distribution, multiplied by the number of cubic nanometers
N_flat=round(33.4*x_length*Radius^2*pi);
N_var=round(.03*randn*N_flat );% 3% variance 
N=N_flat+N_var;
N_H3O=length(x_sulf);
N_H2O = N; 
if mode==2
    tracker = zeros(N_H3O-10,t_max-1);
end

x_H2O=zeros(1,N_H2O);
y_H2O=zeros(1,N_H2O);
z_H2O=zeros(1,N_H2O);


protons=0;
prot_track=zeros(1, t_max);
side_track=zeros(1, t_max);


% Print initial conditions
fprintf('\nRadius: %g nm Length: %g nm \n%g Water molecules, %g Ion pairs, %g cycles \n',Radius, x_length,N, N_H3O,t_max); 



%% Initial placement

ran=1; % 1= random starting placement, 2=specific starting placement

[x_H2O, y_H2O, z_H2O]=coordinate_return3d(x_H2O,y_H2O,z_H2O, x_max, radmax, ran);

order=randperm(N_H2O);
x_temp=x_H2O;
y_temp=y_H2O;
z_temp=z_H2O;
counter=1;
for i=order
   x_H2O(counter)=x_temp(i);
   y_H2O(counter)=y_temp(i);
   z_H2O(counter)=z_temp(i);
   
   counter=counter+1;
end

x_H3O=x_H2O(1:N_H3O);
y_H3O=y_H2O(1:N_H3O);
z_H3O=z_H2O(1:N_H3O);

x_H2O=x_H2O(N_H3O+1:end);
y_H2O=y_H2O(N_H3O+1:end);
z_H2O=z_H2O(N_H3O+1:end);

N_H2O=N-N_H3O;



%% Initialize sulfonates

sulfonate_owner=sulfonate_ownership_function(x_H3O,y_H3O,z_H3O, x_sulf,y_sulf,z_sulf,Radius,x_length,scale,r_trap);

[H2O_grid]=grid_locations(x_H2O,y_H2O,z_H2O,x_H3O,y_H3O,z_H3O,x_max,y_max,z_max);

% sulf_R_grid=5;

%% Calculate E-field at all points
[E_x_grid, sulf_R_grid,sulf_E_grid, P_E_grid,E_grid]...
    =sulf_e_field(x_max,y_max,z_max,x_sulf,y_sulf,z_sulf, radmax,scale, H2O_grid,Radius, E_strength);


%% Initiate loop and Sort molecules


hist_averaging=zeros(t_max,x_max-1); %Axial H3O data
hist_rad=zeros(t_max,radmax-2); % Radial H3O data


%Reorganize H2O molecules into random order.  This helps randomize H3O-H2O
%interactions
order= randperm(N_H2O);%Randomly reorder molecules
x_temp=x_H2O;
y_temp=y_H2O;
z_temp=z_H2O;
counter=1;
for i=order
   x_H2O(counter)=x_temp(i);
   y_H2O(counter)=y_temp(i);
   z_H2O(counter)=z_temp(i);
   
   counter=counter+1;
end
tlapse=toc;

fprintf('%.0f seconds elapsed. Initiate loop\n\n', tlapse)
tic
while t<t_max

    


%% Translation

[H2O_grid]=grid_locations(x_H2O,y_H2O,z_H2O,x_H3O,y_H3O,z_H3O,x_max,y_max,z_max);


[x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, currentd_count,sulf_currentd_count]= ...
    translation3d(x_H2O, y_H2O, z_H2O,x_H3O, y_H3O, z_H3O, x_sulf, y_sulf, z_sulf,...
    x_max, Radius ,sulfonate_owner,H2O_grid, E_grid,r_trap,exp_barrier,mode);


sulfonate_owner=sulfonate_ownership_function(x_H3O,y_H3O,z_H3O, ...
    x_sulf,y_sulf,z_sulf,Radius,x_length,scale,r_trap);



[H2O_grid]=grid_locations(x_H2O,y_H2O,z_H2O,x_H3O,y_H3O,z_H3O,x_max,y_max,z_max);

if mode==2
    for ijk=1:(N_H3O-10)
    tracker(ijk,t) = x_H3O(ijk);
    end
end
if t>1 && mode==2
for ijk=1:N_H3O-10
   if tracker(ijk,t)>tracker(ijk,t-1)+10
       tracker(ijk,t)=tracker(ijk,t-1)-1;
   end
  if tracker(ijk,t)<tracker(ijk,t-1)-10
       tracker(ijk,t)=tracker(ijk,t-1)+1;
   end
end
end
%% Grotthuss Mechanism

[x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O,currentg_count,sulf_currentg_count]=...
    grotthuss3d(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O,x_sulf, y_sulf, z_sulf,... 
    x_max, Radius, P_G,sulfonate_owner,H2O_grid, E_grid, r_trap,exp_barrier);



%% Proton hopping at the borders

[x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, currentb_count,protons,side_track(t)] = ...
    borderprotons3d(x_H2O, y_H2O,z_H2O, x_H3O, y_H3O,z_H3O, x_max,...
    length(x_sulf), target_I,protons, mode);

prot_track(t)=protons;

sulfonate_owner=sulfonate_ownership_function(x_H3O,y_H3O,z_H3O, x_sulf,y_sulf,z_sulf,Radius,x_length,scale,r_trap);


%% Plot 


t=t+1;
if t> t_max+1 && mod(t,5)==0 % Waits until enough timesteps have passed for molecules to diffuse
    print_molecules(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, ...
        x_sulf, y_sulf, z_sulf, ...
        Radius, x_length, scale, printmode,t, sulfonate_owner, r_trap);
end
%% Update current states

%current_counts = counts/timestep
%1 count = e*mark
%current=charge/time; 
%current = current_counts*e*mark/timestep/x_length

mid_currentd(t)=currentd_count/x_length;%Turn counts into current
mid_currentg(t)=currentg_count/x_length;
mid_currentb(t)=currentb_count/x_length;
mid_currents(t)=(sulf_currentd_count+ sulf_currentg_count)/x_length;

if t>1000
    current_d=mean(mid_currentd(t-1000+1:t));
    current_g=mean(mid_currentg(t-1000+1:t));
    current_b=mean(mid_currentb(t-1000+1:t));
    current_s=mean(mid_currents(t-1000+1:t));
    current_tot=current_d+current_g+current_b;
end

if t<1000

    current_d=mean(mid_currentd(1:t));
    current_g=mean(mid_currentg(1:t));
    current_b=mean(mid_currentb(1:t));
    current_s=mean(mid_currents(1:t));
    current_tot=(current_d+current_g+current_b);
end
    
    

%% Arrange data, End timestep

%Radial H3O data
hist_data=zeros(1,radmax-2);
for i=1:length(hist_data)
    hist_data(i)=length(find( sqrt(y_H3O(1,:).^2+z_H3O(1,:).^2)<i+1 & sqrt(y_H3O(1,:).^2+z_H3O(1,:).^2)>=i))*scale^2/(pi*(i+1)^2-pi*(i)^2)/x_length;
end
hist_rad(t,:)=hist_data;


%Axial H3O data
hist_data=zeros(1,x_max-1);
for i=1:length(hist_data)
    hist_data(i)=length(find(x_H3O(1,:)==i))/(pi*Radius^2)*scale/1;
end
hist_averaging(t,:)=hist_data;




%% Approximate completion and time left
if t==1
   disp('First step'); 
end
if mod(t,200)==0   
    fprintf(1,'%5.0f of %.0f   ', t, t_max);
    time_el=toc;
    time_sec=mod(time_el,60);
    time_min_tot=floor(time_el/60);
    time_min=abs(mod(time_min_tot,60));
    time_hour=floor(time_el/3600);
    
    time_left= toc/t * (t_max-t);
    time_l_sec=mod(time_left,60);
    time_l_min_tot=floor(time_left/60);
    time_l_min=mod(time_l_min_tot,60);
    time_l_hour=floor(time_left/3600);
    
    fprintf('Elapsed time: %.0f:%2.0f:%2.0f   ', time_hour,time_min,time_sec);

    fprintf( 'To Go: %.0f:%2.0f:%2.0f  \n',time_l_hour,time_l_min,time_l_sec);
    pause(.1);
 
         
end


end

%% Loop Ends, finish up

rho_plot=sum(hist_averaging)/t; %average concentration

steady_rho=zeros(1,length(rho_plot(1,:))); % Concentration leaving out initialization effects
for i=1:length(steady_rho)
   steady_rho_pre(i)=sum(hist_averaging(end-current_ticks:end,i))/current_ticks; 
end
axial_vector_pre=(1:x_max-1)/scale;

steady_rho=steady_rho_pre(5:end-5);%Hopefully eliminate edge effects
axial_vector=axial_vector_pre(5:end-5);

steady_rho_rad=zeros(1,radmaxmax-2);
for i=1:length(steady_rho_rad)
    if(i>radmax-2)
        steady_rho_rad(i)=0;
    else
        steady_rho_rad(i)=sum(hist_rad(end-3000:end,i))/3001; 
    end
end

grad_rho=gradient(steady_rho);
del_rho=mean(grad_rho);
radial_vector=(1:length(rho_plot))/scale;
steady_radial_vector=(1:length(steady_rho_rad))/scale;

current_d=mean(mid_currentd(end-current_ticks+1:end));
current_g=mean(mid_currentg(end-current_ticks+1:end));
current_b=mean(mid_currentb(end-current_ticks+1:end));
current_s=mean(sulf_currentd(end-current_ticks+1:end))+mean(sulf_currentg(end-current_ticks+1:end));
current_tot=current_d+current_g+current_b;
current_density=current_tot/(pi*Radius^2);


lambda=N/length(x_sulf);
rho_ave=length(x_sulf)/(x_length*pi*Radius^2); %number of ions per nm^3



%% Print results

fprintf('Total current: %g Diffusive current: %g  Grotthuss current: %g Sulf current: %g  \n', current_tot, current_d, current_g, current_s);
%Parameters: sulf_E_field:E_strength
%Grotthuss/Translation: exp_barrier
%P_g
%Radius


% The rest of the program is for viewing and debugging, so it does not run
% with the program normally. It will only run when specifically selected or
% copied.
trackave=zeros(1,t_max-1);

        if mode==2

%%
        for ijk=1:(t_max-1)
            trackave(ijk)=sum(tracker(:,ijk));
        end
        clf;
        figure(1);
        x=(1:length(steady_rho))/scale;
        plot(x, steady_rho);
        xlabel('distance (nm)')
        ylabel('concentration (1/nm^3)')
        pause(.01);



        %% Plot radial concentration

% 
%         figure(2)
%         plot(steady_radial_vector(2:end-3), steady_rho_rad(2:end-3));
%         xlabel('radius (nm)')
%         ylabel('concentration (1/nm^3)')
%         pause(.01);
%         
        figure(3)
        plot(trackave);

        %% Final Plot
% 
%         printmode = 3;
%         print_molecules(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, ...
%                 x_sulf, y_sulf, z_sulf, ...
%                 Radius, x_length, scale, printmode,t, sulfonate_owner,r_trap);        
        end




return;