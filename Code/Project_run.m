disp('New Run-Press Enter to continue');

pause;
clear;


Rad=2;
global t;
t=0;

for i=11:4020   

mode=mod(i,2)+1;

switch mode
    case 2
        E_str=4;
    case 1
        E_str=0;

end
modval=mod(i,21)/20;


[E_g , exp_barrier,  Radius , E_strength,  current_tot  ,...
    current_d ,  current_g , del_rho,  current_ticks, ...
    current_density, lambda, N, rho_ave, steady_rho_rad]=...
    main_diffusion3d_x6(2+modval,E_str, mode);
    

fprintf(1,'\n Completed Run %g \n\n ',i);


srdata=[E_g  (-exp_barrier)  Radius E_strength current_tot ...
    current_d   current_g  del_rho current_ticks current_density...
    lambda N rho_ave];
rhodata=[Radius steady_rho_rad];

if mode==1
    save( 'srdata_MkIII_1.txt','srdata' ,'-ASCII','-append');
end
if mode==2
    save( 'srdata2_MkIII_2.txt','srdata' ,'-ASCII','-append');
end

 save( 'rhodata.txt','rhodata' ,'-ASCII','-append');
     
end

%Play a noise to let me know program has finished
pause(.1);
x=0:3000;
y=sin(.3.*x.*(1/3000.*x));
sound(y);
pause(.01);

