%This program is meant to be used with the radial concentration data from
%the latest versions of the simulationas of August 7 2007


%% Fit model to simulation data

disp('New Fit')

rho_0=.02;
q=.15;
r_off=.15;

y_exp=hist_plot_rad(2:end);
r=2:radmax-1;

y_fit=rho_0./(1-(q*r-r_off).^2/8).^2;
chisq=sum((y_exp-y_fit).^2);


d_rho=.001;
d_q=.01;
d_r=.01;
it=0;
for i=1:14    
    
    chisq_old=chisq+1;
    rho_0_new=rho_0;
    while chisq<chisq_old
        rho_0=rho_0_new;
        chisq_old=chisq;
        rho_0_new=rho_0+d_rho;
        y_fit=rho_0_new./(1-(q*r-r_off).^2/8).^2;
        chisq=sum((y_exp-y_fit).^2);
        it=it+1;
        if mod(it,100)==0
            plot(r,y_fit,r,y_exp)
            pause(.05);
        end
            
    end

    y_fit=rho_0./(1-(q*r-r_off).^2/8).^2;
    chisq=sum((y_exp-y_fit).^2);
    chisq_old=chisq+1;
    q_new=q;
    while chisq<chisq_old
        q=q_new;
        chisq_old=chisq;
        q_new=q+d_q;
        y_fit=rho_0./(1-(q_new*r-r_off).^2/8).^2;
        chisq=sum((y_exp-y_fit).^2);
        it=it+1;
        if mod(it,100)==0
            plot(r,y_fit,r,y_exp)
            pause(.05);
        end
    end

    y_fit=rho_0./(1-(q*r-r_off).^2/8).^2;
    chisq=sum((y_exp-y_fit).^2);
    chisq_old=chisq+1;
    r_new=r_off;
    while chisq<chisq_old
        r_off=r_new;
        chisq_old=chisq;
        r_new=r_off+d_r;
        y_fit=rho_0./(1-(q*r-r_new).^2/8).^2;
        chisq=sum((y_exp-y_fit).^2);
        it=it+1;
        if mod(it,100)==0
            plot(r,y_fit,r,y_exp)
            pause(.05);
        end
    end

    d_rho=-.5*d_rho;
    d_q=-.5*d_q;
    d_r=-.5*d_r;
end


disp('rho_0')
disp(rho_0)
disp('q')
disp(q)
disp('r_off')
disp(r_off)

y_fit=rho_0./(1-(q*r-r_off).^2/8).^2;
plot(r,y_fit,r,y_exp)


%% solve for q, rho_0, based on the model being used

x_length=1;
Radius=1;
scale=10/.5774;
radmax=floor(Radius*scale);
x_max=round(x_length*scale);

kb=1.38*10^-23; %Boltzman's constant, J/K
T=373; %Temperature, K
el=1.602*10^-19; %Electron charge, C
ep_0 = 8.854*10^-12; %Electric constant epsilon_0, C^2/N*m^2 
kappa=55; %Dielectric constant of water
ticks=.0577*10^-9; % Conversion from ticks to meters, m/tick
x_max1=x_length*10^-9;
radmax1=Radius*10^-9;

rho_0=6*10^27; %Concentration at center #/m^3
q_square=el^2*rho_0/(kappa*ep_0*kb*T);
q=sqrt(q_square);
r=(0:98)/100*(radmax1);
x=q*r;
d_rho=.01*10^27;

for i=1:5000
    rho=rho_0./(1-x.^2/8).^2;
    d_r=r(2)-r(1);
    N_tot=sum(2*pi*r.*d_r.*rho*x_max1);
    Vol_tot=(pi*(radmax1)^2*x_max1);
    rho_ave=N_tot/Vol_tot;

    if rho_ave<2*10^27
        rho_0=rho_0+d_rho;
    end
    if rho_ave>2*10^27
        rho_0=rho_0-d_rho;
    end
    if mod(i,1000)==0
        d_rho=1/10*d_rho;
    end

    q_square=el^2*rho_0/(kappa*ep_0*kb*T);
    q=sqrt(q_square);
    r=(0:98)/100*(radmax1);
    x=q*r;
end

plot(r,rho)
disp('rho average')
disp(rho_ave)
disp('end');


