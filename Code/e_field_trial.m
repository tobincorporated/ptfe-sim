


q_enc=zeros(1,length(steady_radial_vector));
for i=2:length(steady_radial_vector)
    q_enc(i)=trapz(steady_radial_vector(1:i),steady_hist_plot_rad(1:i).*steady_radial_vector(1:i)*2*pi*x_length);
end

lambda=q_enc/x_length; % # of charges/nm

N_win=20;
G=gausswin(N_win);
rho_smth_whole=conv(G,steady_hist_plot_rad);
rho_smth=rho_smth_whole(N_win/2:end-N_win/2);
rho_smth=rho_smth*trapz(steady_hist_plot_rad)/trapz(rho_smth);


q_smth=q_enc*0;
for i=2:length(steady_radial_vector)
    q_smth(i)=trapz(steady_radial_vector(1:i),rho_smth(1:i).*steady_radial_vector(1:i)*2*pi*x_length);
end

lam_smth=q_smth/x_length;%ions/nm --> Coulombs/m

R_si=steady_radial_vector;

V_r=0*q_enc;% Radially variant part of the potential
V_E=V_r;%Origianl field used in simulation

for i=1:length(steady_radial_vector)-1
    V_r(end-i)=trapz(R_si(end-i:end),lam_smth(end-i:end)./R_si(end-i:end));
    V_E(i)=E_grid(5,radmax+3,radmax+3+i);
end
V_E(1+i)=E_grid(5,radmax+3,radmax+3+i+1);
V_E=V_E-V_E(1);

%%
T=300;
V_c=1.602*1.602/(2*pi*8.85*1.38*T);
V_exp=-19-19+9+12+23;


V_tot=V_r.*V_c.*10^(V_exp);%*(3*10^-3);
V_tot=V_tot-V_tot(1);

V_P=-abs(2*log(1-(R_si/(Radius)).^2));

%%
plot(R_si,rho_smth);
xlabel('radius (nm)')
ylabel('concentration (1/nm^3)')



%%

semilogy(R_si,rho_smth);
xlabel('radius (nm)')
ylabel('concentration (1/nm^3)')

%%
plot(R_si,steady_hist_plot_rad);
xlabel('radius (nm)')
ylabel('concentration (1/nm^3)')

%%
% 
% V_new=-.15*R_si.^4;
% plot(R_si,V_new,R_si,V_tot)
% xlabel('radius (nm)')
% ylabel('Energy (kb T)')
% 
% V_E=V_new;

%%
plot( R_si,V_E)
xlabel('radius (nm)')
ylabel('Energy (kb T)')

%%
plot(R_si,V_E, R_si,V_tot)
xlabel('radius (nm)')
ylabel('Energy (kb T)')


