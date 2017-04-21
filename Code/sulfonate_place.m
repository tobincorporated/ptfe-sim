function[x_sulf,y_sulf,z_sulf]=sulfonate_place(Radius_si, x_length_si,scale)


%% Initialize
% Radius_si=1;
% x_length_si=10;
scale=10/.5774;  %scale = ticks/nm


%Values from nanometers to ticks
Radius=round(Radius_si*scale);
x_length=round(x_length_si*scale);
circ=round(2*pi*Radius);

%Sulfonate groups are separated by .6-1.1 nm
d_min_si=.6;
d_max_si=1.1;
d_min=round(d_min_si*scale);
d_max=round(d_max_si*scale);

% Flatten out the cylinder wall and turn it into a grid
sulf_grid=zeros(x_length,circ ); 
%3= location of sulfonate, 
%2= too close for another sulfonate, 
%1 = optimum distance for sulfonate, 
%0= completely empty space

[w_0,l_0]=find(sulf_grid==0);
[w_1,l_1]=find(sulf_grid==1);

o_filled = 0;
empty=1;
%% determine new sulfonate

while o_filled==0 % As long as there is unoptimized space on the cylinder
    
    switch empty

        case 0 %Pick a random spot where the grid is 1 to place a new sulfonate
            spots=length(w_1); 
            index=ceil(spots*rand);
            w_n=w_1(index);
            l_n=l_1(index);
            
            
        case 1 %Pick a random spot where the grid is 0 to place a new sulfonate
            spots=length(w_0);
            index=ceil(spots*rand);
            w_n=w_0(index);
            l_n=l_0(index);


    end
    
    
%% Set surroundings of sulfonates

    sulf_grid(w_n,l_n)=3;
    
    %Set optimum sulfonate locations
    w_min=w_n-d_max;
    if w_min<1
        w_min=1;
    end    
    w_max=w_n+d_max;
    if w_max>x_length
        w_max=x_length;
    end    
    l_min=l_n-d_max;
    if l_min<1
        l_min=1;
    end
    l_max=l_n+d_max;
    if l_max>circ
        l_max=circ;
    end

    
    for i=w_min:1:w_max
        for j=l_min:1:l_max
            
            d_i=sqrt( (w_n-i)^2 + (l_n-j)^2);
            
            %Set nearby spots to 1
            if d_i < d_max && sulf_grid(i,j )< 1
                sulf_grid(i,j)=1;
            end
            
            %Set closer spots to 2
            if d_i < d_min && sulf_grid(i,j) < 2
                sulf_grid(i,j)=2;
            end
        end
    end
         
%% wrap-around

    %Set surroundings for wrap-around, since circumference wraps at
    %endpoints
    
    if abs(circ-l_n)<d_max %If l_n is close to the wrap edge
        wrap=1; %set whether to wrap up or down
        l_min=circ-d_max;
        l_max=circ;
        if l_n>circ/2
            wrap=-1;
            l_min=1;
            l_max=d_max;
        end

        l_n=l_n+wrap*circ;
        for i=w_min:1:w_max
            for j=l_min:1:l_max

                d_i=sqrt( (w_n-i)^2 + (l_n-j)^2);

                if d_i < d_max && sulf_grid(i,j )< 1
                    sulf_grid(i,j)=1;
                end
                if d_i < d_min && sulf_grid(i,j) < 2
                    sulf_grid(i,j)=2;
                end
            end
        end
    end
    
    
    
   
    
    [w_3,l_3]=find(sulf_grid==3);%Array of actual sulfonate locations on grid
    [w_2,l_2]=find(sulf_grid==2);
    [w_1,l_1]=find(sulf_grid==1);
    [w_0,l_0]=find(sulf_grid==0);
    o_filled = isempty(w_1);
    empty=0;
    
%% Print current map
%     clf
%     hold on;
%     plot(w_3/scale,l_3/scale,'.k','MarkerSize',24)
%     plot(w_2/scale,l_2/scale,'+r')
%     plot(w_1/scale,l_1/scale,'+b')
%     plot(w_0/scale,l_0/scale,'+g')
%     axis equal
%     hold off;
%     pause();

%% Set values 

    x_sulf=w_3'; 
    theta_sulf=(l_3')./(2*pi);
    y_sulf=round(Radius*sin(theta_sulf));
    z_sulf=round(Radius*cos(theta_sulf));



end


%% Test grid
% Test_grid=0*sulf_grid;
% ions=length(w_3);
% for i=1:x_length
%     for j=1:circ    
%         for i2=1:ions           
%             Test_grid(i,j)= Test_grid(i,j)+1/((w_3(i2)-i)^2+(l_3(i2)-j)^2);
%             if w_3(i2)==i && l_3(i2)==j
%                 Test_grid(i,j)=1.5;
%             end
%         end
%     end
% end
% lnTest= log(Test_grid);
% %contour(Test_grid,30);
% contour(lnTest,30);
% axis equal
return;