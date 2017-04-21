function [x_H2O, y_H2O, z_H2O] = coordinate_return3d(x_H2O,y_H2O,z_H2O, x_max, radmax, ran)




N_H2O=length(x_H2O);

switch ran
    case 1
        
        %Create a smaller matrix with fewer elements that are each spaced
        %out further, making the molecules spaced well and avoiding
        %accidental overlap
        x_cut=floor(x_max/3);
        y_cut=floor(2*radmax/3);% *2 because rad==> diameter
        z_cut=floor(2*radmax/3);    
        H2O_grid=zeros(x_cut,y_cut,z_cut); %3-d grid of 0's. a 1 signifies the location of a H2O
        
        %Following loop seeds the grid with N H2O centers
        while sum(sum(sum(H2O_grid))) < N_H2O % there should be as many 1's as there are H2O's
            y_rand=(ceil(rand*y_cut));
            z_rand=(ceil(rand*z_cut));
            if (y_rand-y_cut/2)^2+(z_rand-z_cut/2)^2 < radmax^2/9-2
                x_rand=(ceil(rand*x_cut));
                H2O_grid(x_rand,y_rand,z_rand)=1;
            end
        end
        
        points=find(H2O_grid); %returns a 1-D array of points that has to be reinterpreted for x,y,z values
        x= (mod(points, x_cut))*3+1;
        y=(ceil ( mod(points,x_cut*y_cut) /x_cut)-round(y_cut/2))*3;
        z=(ceil(points./y_cut./x_cut)-round(z_cut/2))*3;
        
        %Oxygen location
        x_H2O(:)=x;
        y_H2O(:)=y;
        z_H2O(:)=z;
        
       
        
        
    case 2
        % The following nested loops place water molecules in the x-center plane,
        % creating a full circular cross-section of molecules.  This creates a
        % consistent, simple starting point
        index=1;
        y_min=-radmax;
        y_max=radmax;
        z_min=-radmax;
        z_max=radmax;
        for i=round(x_max/2)-1:1:x_max-1
            for j=y_min+2:3:y_max-1
                for k=z_min+2:3:z_max-1
                    if index<=N_H2O && k^2+j^2 <= (radmax-2)^2
                           x_H2O(index)=i;
                           y_H2O(index)=j;
                           z_H2O(index)=k;


                           index=index+1;

                    end
                end
            end
        end
end



% %Test print the created molecules
% x_O2=x_H2O(1,:)/scale;
% y_O2=y_H2O(1,:)/scale;
% z_O2=z_H2O(1,:)/scale;
% 
% x_H2= [x_H2O(2,:) x_H2O(3,:)]/scale;
% y_H2= [y_H2O(2,:) y_H2O(3,:)]/scale;
% z_H2= [z_H2O(2,:) z_H2O(3,:)]/scale;
% 
% x_H2Oline=zeros(3,length(x_H2O(1,:)));
% x_H2Oline(1,:)=x_H2O(2,:)/scale;
% x_H2Oline(2,:)=x_H2O(1,:)/scale;
% x_H2Oline(3,:)=x_H2O(3,:)/scale;
% y_H2Oline=zeros(3,length(x_H2O(1,:)));
% y_H2Oline(1,:)=y_H2O(2,:)/scale;
% y_H2Oline(2,:)=y_H2O(1,:)/scale;
% y_H2Oline(3,:)=y_H2O(3,:)/scale;
% z_H2Oline=zeros(3,length(x_H2O(1,:)));
% z_H2Oline(1,:)=z_H2O(2,:)/scale;
% z_H2Oline(2,:)=z_H2O(1,:)/scale;
% z_H2Oline(3,:)=z_H2O(3,:)/scale;
% 
% clf;
% hold on
% plot3(x_O2, y_O2, z_O2,'bo', 'MarkerSize',12);
% plot3(x_H2,y_H2,z_H2,'b.', 'MarkerSize',24) ;
% hold off
% line(x_H2Oline,y_H2Oline, z_H2Oline,'Color','b');
% xlabel('x-axis');
% ylabel('y-axis');
% zlabel('z-axis');
% axis equal
% view(90,0);



return;
