Kap=1;
ep_0=8.85*10^-12;
el=1.6*10^-19;
kb=1.38*10^-23;
T=300;
rho_H=33.4*10^27;
R=(2:.01:3) *10^-9;
r=(0:.01:3)*10^-9;
m=22/3*10^9;
rho_0=8*ep_0*Kap*rho_H/m./R./(8*ep_0*Kap+1/kb/T*el^2.*R*rho_H/m);

ep=8*ep_0*Kap;
gamma=el^2*rho_H/(kb*T*m);
rho=zeros(length(R),length(r));

for i=1:length(R)
    for j=1:length(r)
        rho_num=ep*rho_H/m/R(i);
        rho_dem1=gamma*r(j)^2/R(i);
        rho_dem2=ep+gamma*R(i);
        rho(i,j)=rho_num/(1-rho_dem1/rho_dem2)^2;

%         rho(i,j)=5;

    end
end

rho_ave=rho_H./(m.*R);
qsq=-8/kb/T*el^2.*rho_ave./(8*Kap*ep_0-el^2.*R.^2/kb^2/T^2.*rho_ave);

r=0:.01:3;
R=2:.01:3;

plot_rho=log(rho);

for i= 1:length(R)
    for j=1:length(r)
        if r(j)>R(i)
            plot_rho(i,j)=min(min(plot_rho));
%             plot_rho(i,j)=0;
        end
    end
end


surf(r,R,plot_rho, plot_rho  );
xlabel('radius r(nm)');
ylabel('radius R(nm)');
axis equal;
axis([0 3 2 3]);
colorbar ;
view(90,90);
shading interp;