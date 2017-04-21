%%
%Zachary Tobin
%Mentor: Philip Taylor
%Case Western Reserve University
%Department of Physics
%SURES program for study of Proton Diffusion in hydrated Nafion membrane

%This program simplifies proton diffusion to a 2-dimensional square
%   lattice.
%In its current state, it will ignore most intermolecular interactions
%   and focus primarily on random diffusive movement, orientation, and
%   Grotthus mechanism
%It is my intention to expand this simulation into 3-dimensions and
%   eventually extend the scope to incorporate more physical forces to be a
%   more realistic model
clear;

for q=2:2 
for force_length=4:4
t_max=2000;

frac_H2O=1;
conc_H3O= .15;
x_length=50;
y_length=50;
P_T   =  1;   %probability during a given timestep that a water molecule will translate
P_R   = 1;  %probability during a given timestep that a water molecule will rotate
P_G   = .8; %probability during a given timestep that a water molecule adjacent to a
             %   hydronium hydrogen will undergo Grotthus hopping

             
             
%% This cell is for initial placement of molecules
            

N_H2O=(100*frac_H2O+round(0*randn)) *x_length/50 * y_length/50;

% N_H3O=round(N_H2O*conc_H3O);
N_H3O=30 *x_length/50  *y_length/50;


x_H2Omas=zeros(5,N_H2O);
y_H2Omas=zeros(5,N_H2O);
x_H3Omas=zeros(6,N_H3O);
y_H3Omas=zeros(6,N_H3O);


for i=1:N_H2O
    proceed=0;
    while proceed==0
        x_H2Omas(1,i)=floor(x_length*rand/3)*3+1;
        x_H2Omas(3,i)=x_H2Omas(1,i);
        x_H2Omas(5,i)=x_H2Omas(1,i);
        y_H2Omas(1,i)=floor(y_length*rand/3)*3+1;
        y_H2Omas(3,i)=y_H2Omas(1,i);
        y_H2Omas(5,i)=y_H2Omas(1,i);
        proceed=1;
        if i>1
            for j=1:i-1
                if x_H2Omas(1,i)==x_H2Omas(1,j) && y_H2Omas(1,i)==y_H2Omas(1,j)
                    proceed=0;
                end
            end
        end
    end    
end




for i=1:N_H3O
    proceed=0;
    while proceed==0
        x_H3Omas(1,i)=floor(x_length*rand/3)*3+1;
        x_H3Omas(3,i)=x_H3Omas(1,i);
        x_H3Omas(5,i)=x_H3Omas(1,i);
        y_H3Omas(1,i)=floor(y_length*rand/3)*3+1;
        y_H3Omas(3,i)=y_H3Omas(1,i);
        y_H3Omas(5,i)=y_H3Omas(1,i); 
        proceed=1;
        if i>1
            for j=1:i-1
                if x_H3Omas(1,i)==x_H3Omas(1,j) && y_H3Omas(1,i)==y_H3Omas(1,j)
                    proceed=0;
                end 
                d= sqrt( (x_H3Omas(1,i)- x_H3Omas(1,j))^2 + (y_H3Omas(1,i)- y_H3Omas(1,j))^2);
                if d<force_length
                    proceed=0;
                end
            end
            
        end
        
        for j=1:N_H2O
            if x_H3Omas(1,i)==x_H2Omas(1,j) && y_H3Omas(1,i)==y_H2Omas(1,j)
                proceed=0;
            end
        end 
        
    end
end

x_H2O=x_H2Omas(1,:);
y_H2O=y_H2Omas(1,:);
x_H3O=x_H3Omas(1,:);
y_H3O=y_H3Omas(1,:);




[x_H2Omas, y_H2Omas] = coordinate_return2d(x_H2Omas, y_H2Omas, 2);
[x_H3Omas, y_H3Omas] = coordinate_return2d(x_H3Omas, y_H3Omas, 3);




x_H2=[ x_H2Omas(2,:) x_H2Omas(4,:)] ;
y_H2=[ y_H2Omas(2,:) y_H2Omas(4,:)] ;

x_H3=[x_H3Omas(2,:) x_H3Omas(4,:) x_H3Omas(6,:)] ;
y_H3=[y_H3Omas(2,:) y_H3Omas(4,:) y_H3Omas(6,:)] ;



% clf;
% subplot(1,2,1);
% plot(x_H2O,y_H2O,'bo',x_H3O,y_H3O,'ro',x_H2,y_H2,'b.',x_H3,y_H3,'r.');
% line(x_H2Omas,y_H2Omas,'Color','b');
% line(x_H3Omas,y_H3Omas,'Color','r');
% axis([ 0 x_length 0 y_length]);
% 
% subplot(1,2,2);         
hist_data=zeros(1,x_length-9);
for i=1:length(hist_data)
    hist_data(i)=length(find(x_H3Omas(1,:)>=i & x_H3Omas(1,:)<=9+i));
end
% plot(5:x_length-5, hist_data);
% axis([0 x_length 0 length(x_H3O)/2]);


ave_new=0;
ave_old=100;
equil=0;
hist_averaging=zeros(2000,x_length-9);
tic
for t=1:t_max


%% This cell is for Random Translation in the plane
H2O_master(1:5,:)=x_H2Omas(:,:);
H2O_master(6:10,:)=y_H2Omas(:,:);
H2O_master=sortrows(H2O_master')' ;
x_H2Omas(:,:)=H2O_master(1:5,:);
y_H2Omas(:,:)=H2O_master(6:10,:);

H3O_master(1:6,:)=x_H3Omas(:,:);
H3O_master(7:12,:)=y_H3Omas(:,:);
H3O_master=sortrows(H3O_master')';
x_H3Omas(:,:)=H3O_master(1:6,:);
y_H3Omas(:,:)=H3O_master(7:12,:);

translation_coord = translation_backup(P_T, x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas, x_length, y_length, 2, force_length) ;
for i=1:length(translation_coord(1,:,1))
    for j=1:length(translation_coord(1,1,:))
        x_H2Omas(i,j)=translation_coord(1,i,j);
        y_H2Omas(i,j)=translation_coord(2,i,j);
    end
end
translation_coord = translation_backup(P_T, x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas, x_length, y_length, 3, force_length) ;
for i=1:length(translation_coord(1,:,1))
    for j=1:length(translation_coord(1,1,:))
        x_H3Omas(i,j)=translation_coord(1,i,j);
        y_H3Omas(i,j)=translation_coord(2,i,j);
    end
end
 





%% This cell is for Random Rotation in the plane



rotate_coord = rotation_backup(P_R, x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas, x_length, y_length, 2);
for i=1:length(rotate_coord(1,:,1))
    for j=1:length(rotate_coord(1,1,:))
        x_H2Omas(i,j)=rotate_coord(1,i,j);
        y_H2Omas(i,j)=rotate_coord(2,i,j);
    end
end
rotate_coord = rotation_backup(P_R, x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas,x_length, y_length, 3);
for i=1:length(rotate_coord(1,:,1))
    for j=1:length(rotate_coord(1,1,:))
        x_H3Omas(i,j)=rotate_coord(1,i,j);
        y_H3Omas(i,j)=rotate_coord(2,i,j);
    end
end





%% This cell is for the Grotthuss mechanism



[x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas] = ...
    grotthus2d(P_G, x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas);





%% For maintaining left concentration

% [x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas] = ...
%     fixedside2d(x_H2Omas, y_H2Omas, x_H3Omas, y_H3Omas, 10, 10*q );
% 




%% This cell is for printing out the results of the timestep's movement
x_H2O=x_H2Omas(1,:);
y_H2O=y_H2Omas(1,:);
x_H2=[ x_H2Omas(2,:) x_H2Omas(4,:)] ;
y_H2=[ y_H2Omas(2,:) y_H2Omas(4,:)] ;
x_H3O=x_H3Omas(1,:);
y_H3O=y_H3Omas(1,:);
x_H3=[ x_H3Omas(2,:) x_H3Omas(4,:) x_H3Omas(6,:)] ;
y_H3=[ y_H3Omas(2,:) y_H3Omas(4,:) y_H3Omas(6,:)] ;

% pause(.2);
% % subplot(1,2,1);
% plot(x_H2O,y_H2O,'bo',x_H3O,y_H3O,'ro',x_H2,y_H2,'b.',x_H3,y_H3,'r.');
% line(x_H2Omas,y_H2Omas,'Color','b');
% line(x_H3Omas,y_H3Omas,'Color','r');
% axis([ 0 x_length 0 y_length]);


hist_data=zeros(1,x_length-9);
for i=1:length(hist_data)
    hist_data(i)=length(find(x_H3Omas(1,:)>=i & x_H3Omas(1,:)<=9+i));
end
hist_averaging(t,:)=hist_data;
% hist_plot=sum(hist_averaging);
% subplot(1,2,2);
% plot(5:x_length-5, hist_plot/t);
% axis([ 0 x_length 0 20]);
if mod(t,1000)==0
    toc
end



end

hist_plot=sum(hist_averaging)/t;
subplot(1,2,2);
plot(5:x_length-5, hist_plot);
axis([ 0 x_length 0 20]);

time=num2str(t);
width=num2str(x_length);
height=num2str(y_length);
f_length=num2str(force_length);
conc_set=num2str(5*q);
savename=[height 'by' width 'times' time 'with_length' f_length 'and_conc' conc_set];
save(savename, 'hist_plot');


end
end