%%

% clf;
colormap('default');

vi= 2; %1 = normal view, 2= top view
grid_num =4 ;
switch grid_num
    case 1
        grid_use=P_E_grid;
    case 2
        grid_use=sulf_E_grid;
    case 3
        grid_use=sulf_pos_E_grid;
    case 4
        grid_use=-E_grid;
%             grid_use=sulf_E_grid+P_E_grid+sulf_pos_E_grid;
    case 5
        grid_use=zeros(length(E_grid(1,:,1)));
        for i=1:length(grid_use)
            for j=1:length(grid_use)
               grid_use(i,j)=-E_grid(5,i,j);
            end
        end
end


height=max(max(max(abs(grid_use))));
for t=0:0
    if t==radmax;
        
    end
    radius=radmax-t;
    dx=1;
    x=round(1:dx:x_max+6);
    dtheta=dx/radius;
    if dtheta>2*pi
        dtheta=2*pi;
    end
    theta=0:dtheta:2*pi;


    y=round(radius*sin(theta))+radmax+3;
    z=round(radius*cos(theta))+radmax+3;
    newgrid=zeros(x_max+6,length(theta));
    for i=1:length(x)
        for j=1:length(theta)
            newgrid(i,j)=abs(grid_use(x(i),y(j),z(j)));
            if newgrid(i,j) ==inf
               newgrid(i,j)=4; 
            end
        end
    end
%     newgrid=log(newgrid);

%     newgrid=newgrid*32/height; %*x_max/height;
    
%     figure(1);
%     surf(theta*radius/scale,x/scale,newgrid/scale)
%     axis([0 2*pi*radius 0 x_max 0 height])
    surf(theta,x/scale,newgrid/scale)
    axis([0 2*pi 0 x_max/scale 0 height/scale])
    
    xlabel(texlabel('theta'));
    ylabel('axis length');
    colorbar
    title(radius/scale);
%     axis equal

switch vi
    case 1
        view(40,20)
    case 2
        view(0,90)
end

    x_val=x_sulf+3;
    th_val=y_sulf+radmax+3;
    for i=1:length(y_sulf)
        t_th0=atan(abs(y_sulf(i)/z_sulf(i)));
        if y_sulf(i)>=0 && z_sulf(i)>=0
            th_val(i)= t_th0;
        end
        if y_sulf(i)>=0 && z_sulf(i)<=0
            th_val(i)=pi- t_th0;
        end
        if y_sulf(i)<=0 && z_sulf(i)<=0
            th_val(i)= pi+t_th0;
        end
        if y_sulf(i)<=0 && z_sulf(i)>=0
            th_val(i)=2*pi- t_th0;
        end
        if th_val(i)>2*pi
            pause;
        end

    end



    time=num2str(t);
    saveas(gcf, ['animation/animation' time],'jpg');
end
% x_val=x_val/scale;
% th_val=th_val;
% figure(2)
% scatter(th_val,x_val)
% axis([0 2*pi 0 (x_max+3)/scale])


%%
if t>t_max+300
%%

        colormap('default')
        grid_use=zeros(length(E_coarse(1,:,1)));
        new_grid=E_grid;
%         new_grid=smooth3(E_coarse,'gaussian',5,5);
%         new_grid=smooth3(E_grid,'box',3);
        
        
        for i=1:length(grid_use)
            for j=1:length(grid_use)
               grid_use(i,j)=-new_grid(10,i,j);
            end
        end
        pause(.1)
        surf((grid_use))
        pause(.1)
        surf((grid_use))
        pause(.1)
        view(45,0)
        pause(.1)
        axis([0 80 0 80 0 10])
        pause(.1)
        
%%
end