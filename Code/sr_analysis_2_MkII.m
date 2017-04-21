Data1 = load('srdata_MkIII_1.txt');
Data2 = load('srdata2_MkIII_2.txt'); 
Data =  Data1;
En_G=Data(:,1);
En_S=Data(:,2);
Radius=Data(:,3);
En_x=Data(:,4);
I_tot=Data(:,5);
I_dif=Data(:,6);
I_grot=Data(:,7);
grad_rho=Data(:,8);
counts=Data(:,9);
current_dens=Data(:,10);

datapoints=length(En_G);
%% Choose analytic variable


avar_x=3;
avar_y=11;
%     case 1 Grotthuss Energy
%     case 2 Sulfonate Well
%     case 3 Radius/hydration
%     case 4 E_x strength
%     case 5 Total current
%     case 6 diffusive current
%     case 7 Grotthuss Current
%     case 8 del rho
%       9 ticks
%       10 current density
%       11 lambda hydration
%       12 # of h2o molecules
%       13 rho average of H+ ions

    
Data = sortrows(Data, avar_x);

%% Use Straight data

x_all=Data(:,avar_x);
y_all=Data(:,avar_y);



x_data=x_all;y_data=y_all;

%% Average  y measurements at each x
num_x=0;%number of different x-values
now_x=0; % current x-value

%Determine where data values change
for i=1:length(x_data)
    if x_data(i)~=now_x
        num_x=num_x+1;
        now_x=x_data(i);
    end
end
x_mark=ones(1,num_x);% Keeps track of where x-values change
x_ave=zeros(1,num_x);
y_ave=zeros(1,num_x);
y_err=zeros(1,num_x);
x_ind=1;

now_x=x_data(1);% current x-value
for i=2:length(x_data)
    if x_data(i)~=now_x
        x_mark(x_ind)=i-1;
        now_x=x_data(i);
        x_ind=x_ind+1;
    end
end
x_mark(end)=length(x_data);

dev_temp=0;
num_dev=0;
x_ind=1;

%Average values at each point
for i=2:(length(x_mark))
    x_ave(x_ind)=x_data(x_mark(i-1));
    y_ave(x_ind)=mean(y_data(x_mark(i-1)+1:x_mark(i)));
    dev_temp=std(y_data(x_mark(i-1)+1:x_mark(i)));
    num_dev=length(y_data(x_mark(i-1)+1:x_mark(i)));
    y_err(x_ind)=dev_temp/sqrt(num_dev);
    x_ind=x_ind+1;
end
x_ave(end)=x_data(x_mark(end));
errorbar(x_ave(1:end-1),y_ave(1:end-1),y_err(1:end-1));
disp('Analysis Run');