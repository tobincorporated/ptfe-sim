function[]= print_molecules(x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, ...
     x_sulf, y_sulf, z_sulf, ...
    Radius, x_length, scale, printmode,t, sulfonate_owner,r_trap)

%All parameters as they are defined in the main program


%% Plot setup

%% H2O
%Location of O in H2O
x_O2=x_H2O(1,:)/scale;
y_O2=y_H2O(1,:)/scale;
z_O2=z_H2O(1,:)/scale;



%% H3O
%Location of O in H3O
x_O3=x_H3O(1,:)/scale;
y_O3=y_H3O(1,:)/scale;
z_O3=z_H3O(1,:)/scale;

%% Sulfonates

%Sulfonate bonds
x_ownership=zeros(2,length(sulfonate_owner(1,:)));
x_i_sulf=x_sulf(sulfonate_owner(2,:));
x_i_h3o=x_O3(sulfonate_owner(1,:));
x_ownership(1,:)=x_i_h3o;
x_ownership(2,:)=x_i_sulf./scale;

y_ownership=zeros(2,length(sulfonate_owner(1,:)));
y_i_sulf=y_sulf(sulfonate_owner(2,:));
y_i_h3o=y_O3(sulfonate_owner(1,:));
y_ownership(1,:)=y_i_h3o;
y_ownership(2,:)=y_i_sulf./scale;

z_ownership=zeros(2,length(sulfonate_owner(1,:)));
z_i_sulf=z_sulf(sulfonate_owner(2,:));
z_i_h3o=z_O3(sulfonate_owner(1,:));
z_ownership(1,:)=z_i_h3o;
z_ownership(2,:)=z_i_sulf./scale;

d_sulf_2=(x_ownership(1,:)-x_ownership(2,:)).^2+...
    (y_ownership(1,:)-y_ownership(2,:)).^2+...
    (z_ownership(1,:)-z_ownership(2,:)).^2;

d_sulf=sqrt(d_sulf_2);

n_sulf_bond=find(d_sulf<r_trap);
sulf_ind=0;
x_sulf_bond=zeros(2,length(n_sulf_bond));
y_sulf_bond=zeros(2,length(n_sulf_bond));
z_sulf_bond=zeros(2,length(n_sulf_bond));

for i=n_sulf_bond
    sulf_ind=sulf_ind+1;
    x_sulf_bond(:,sulf_ind)=x_ownership(:,i);
    y_sulf_bond(:,sulf_ind)=y_ownership(:,i);
    z_sulf_bond(:,sulf_ind)=z_ownership(:,i);
    
end

%Sulfonate locations
x_sulf_print=x_sulf./scale;
y_sulf_print=y_sulf./scale;
z_sulf_print=z_sulf./scale;


%% Start of the print sections

%It is best to only use one of the print sections at a time. The first is
%best for debugging the program, the second is better for presentation

%% Print

switch printmode
    case 1
        
        clf;
        hold on;
        plot3(x_O2, y_O2, z_O2,'bo', 'MarkerSize',12);
        plot3(x_O3, y_O3, z_O3,'ro', 'MarkerSize',12);
    	line(x_ownership,y_ownership,z_ownership,'Color','m')
        plot3(x_sulf_print,y_sulf_print,z_sulf_print,'.m');


        [z,y,x]=cylinder(Radius+1/scale,31);
        x(2,:)=x_length;
        y=[y(:,end*3/4:end) y(:,1:end*3/4-1)];
        z=[z(:,end*3/4:end) z(:,1:end*3/4-1)];
        surf(x(:,10:end),y(:,10:end),z(:,10:end),z(:,10:end).*0+1)
        colormap([ 0 1 0]);

        
        % Make plot more presentable
        
        shading interp; 
        lighting gouraud; 
        camlight infinite;
        hold off;
        
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        axis([0 x_length -Radius Radius -Radius Radius]);
        axis equal
        view(20,30);
        pause(.2);

%% Print without H2O

case 2
        
        clf;
        hold on;
        plot3(x_O3, y_O3, z_O3,'ro', 'MarkerSize',12);
    	line(x_ownership,y_ownership,z_ownership,'Color','m')
        plot3(x_sulf_print,y_sulf_print,z_sulf_print,'.m');

        hold off;
        
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        axis([0 x_length -Radius Radius -Radius Radius]);
        axis equal
        view(20,30);
        pause(.2);


%% Print 3d for save--very processor-intensive

    case 3

        
        
        
        pause(.05);
        clf('reset');
        hold on;

        % Plot atoms as spheres

        n=3;%precision
        [x, y, z]=sphere(n); 


        %H2O oxygen
        %%radius of oxygen is actually ~60 pm
        r=.04;
        for i=1:length(x_O2);
            x_plot=x_O2(i)+x*r;
            y_plot=y_O2(i)+y*r;
            z_plot=z_O2(i)+z*r;
%            surf(x_plot, y_plot, z_plot, z_plot.*0+1)
            surf(x_plot, y_plot, z_plot, z_plot.*0+1,'FaceAlpha',.2);
        end


        n=6;%precision
        [x, y, z]=sphere(n);      
        
        
        %H3O oxygen
        r=.04;
        for i=1:length(x_O3);
            x_plot=x_O3(i)+x*r;
            y_plot=y_O3(i)+y*r;
            z_plot=z_O3(i)+z*r;
            surf(x_plot, y_plot, z_plot, z_plot.*0+2);
        end

        % Sulfonates
        r=.045;
        for i=1:length(x_sulf);
            x_plot=x_sulf_print(i)+x*r;
            y_plot=y_sulf_print(i)+y*r;
            z_plot=z_sulf_print(i)+z*r;
            surf(x_plot, y_plot, z_plot,  z_plot.*0+3);
        end

        %cylinder
        [z,y,x]=cylinder(Radius+1/scale,31);
        x(2,:)=x_length;
        y=[y(:,end*3/4:end) y(:,1:end*3/4-1)];
        z=[z(:,end*3/4:end) z(:,1:end*3/4-1)];
        
       surf(x,y,z,z.*0+4,'FaceAlpha','interp',...
        'AlphaDataMapping','scaled',...
        'AlphaData',abs((y+Radius)));

         %Cap/proton wall
         
         theta=0:.2:2*pi;
         x=0*theta;
         y=(Radius+1/scale)*sin(theta);
         z=(Radius+1/scale)*cos(theta);

         fill3(x,y,z,z.*0+5,'FaceAlpha',.1)
                
        % Make plotted spheres more presentable
        colormap([0 0 1; 1 0 0 ;1 0 0 ; 1 1 0 ; 1 1 0 ; 0 1 0; 1 0 1]);
        shading interp; 
        lighting gouraud; 
        camlight infinite;

        %plot non-sphere objects
        line(x_sulf_bond,y_sulf_bond,z_sulf_bond,'Color',[.5 0 .5], 'LineWidth', 1);
        hold off;
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        axis equal;
        axis([0 x_length -Radius-.06 Radius+.06 -Radius-.06 Radius+.06]);
        axis off

        view(20,20);
        zoom(1.5)

        time=num2str(t);
        saveas(gcf, ['animation/animation' time],'tif');
        pause(1);


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        
        
        
        
%% Print 3d without H2O

 case 4

        
        
        pause(.05);
        clf('reset');
        hold on;

        % Plot atoms as spheres

        n=7;%precision
        [x, y, z]=sphere(n);      
        
        
        %H3O oxygen
        r=.04;
        for i=1:length(x_O3);
            x_plot=x_O3(i)+x*r;
            y_plot=y_O3(i)+y*r;
            z_plot=z_O3(i)+z*r;
            surf(x_plot, y_plot, z_plot, z_plot.*0+2);
        end


        % Sulfonates
        r=.045;
        for i=1:length(x_sulf);
            x_plot=x_sulf_print(i)+x*r;
            y_plot=y_sulf_print(i)+y*r;
            z_plot=z_sulf_print(i)+z*r;
            surf(x_plot, y_plot, z_plot,  z_plot.*0+3);
        end


        [z,y,x]=cylinder(Radius+1/scale,31);
        x(2,:)=x_length;
        y=[y(:,end*3/4:end) y(:,1:end*3/4-1)];
        z=[z(:,end*3/4:end) z(:,1:end*3/4-1)];
        
       surf(x,y,z,z.*0+4,'FaceAlpha','interp',...
        'AlphaDataMapping','scaled',...
        'AlphaData',abs((y+Radius)));

    
                 %Cap/proton wall
         
         theta=0:.2:2*pi;
         x=0*theta;
         y=(Radius+1/scale)*sin(theta);
         z=(Radius+1/scale)*cos(theta);

         fill3(x,y,z,z.*0+5,'FaceAlpha',.1)
                
                
                
        % Make plotted spheres more presentable
        colormap([1 0 0;  0 1 0 ; 1 1 0; 0 1 0 ; 0 1 0; 0 1 0; 1 0 1]);
        shading interp; 
        lighting gouraud; 
        camlight infinite;

        %plot non-sphere objects
        line(x_sulf_bond,y_sulf_bond,z_sulf_bond,'Color',[.5 0 .5], 'LineWidth', 1);

        hold off;
        xlabel('x-axis');
        ylabel('y-axis');
        zlabel('z-axis');
        axis equal;
        axis([0 x_length -Radius-.06 Radius+.06 -Radius-.06 Radius+.06]);
        axis off

        view(20,20);
        pause(1)
        zoom(1.5)

        time=num2str(t);
        saveas(gcf, ['animation/animation' time],'tif');
        pause(1);


end



return;
