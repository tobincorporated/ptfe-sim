%Data1 = load('Current_Data.txt'); % Original data
Data2 = load('srdata_MkII.txt'); %Newer data
Data =  Data2;
En_G=Data(:,1);
En_S=Data(:,2);
Radius=Data(:,3);
En_x=Data(:,4);
I_tot=Data(:,5);
I_dif=Data(:,6);
I_grot=Data(:,7);
grad_rho=Data(:,8);
counts=Data(:,9);
I_norm=I_tot./counts;
I_normd=I_dif./counts;
I_normg=I_grot./counts;
current_dens=Data(:,10);
Data(:,11)=I_norm;
Data(:,12)=I_normg;
Data(:,13)=I_normd;
datapoints=length(En_G);
%% Choose analytic variable

%Data1(:,9)=Data1(:,5)./(Data1(:,8));

avar_x=3;
avar_y=9;
%     case 1 Grotthuss Energy
%     case 2 Sulfonate Well
%     case 3 Radius/hydration
%     case 4 E_x strength
%     case 5 Total current
%     case 6 diffusive current
%     case 7 Grotthuss Current
%     case 9 Norm current

    
Data = sortrows(Data, avar_x);

%% Use Straight data

x_all=Data(:,avar_x);
y_all=Data(:,avar_y);




%% Determine viable data points

E_str=0.5;
viable_det=ones(1,datapoints);
for i=1:viable_det
    if En_x(i)==E_str
        viable_det(i)=1;
    else viable_det(i)=sqrt(-1);
    end
end

viable=find(viable_det==1);
n_vi=sum(real(viable_det));
x_data=ones(1,n_vi);
y_data=ones(1,n_vi);

for i=1:n_vi
   x_data(i)= x_all(viable(i));
   y_data(i)= y_all(viable(i));
end


%% Average  y measurements at each x
num_x=0;
now_x=0;

%Determine where data values change
for i=1:length(x_data)
    if x_data(i)~=now_x
        num_x=num_x+1;
        now_x=x_data(i);
    end
end
x_mark=ones(1,num_x+1);
x_ave=zeros(1,num_x);
y_ave=zeros(1,num_x);
y_err=zeros(1,num_x);

now_x=x_data(1);
x_ind=2;
for i=2:length(x_data)
    if x_data(i)~=now_x
        x_mark(x_ind)=i-1;
        x_ind=x_ind+1;
        now_x=x_data(i);
    end
end
x_mark(end)=length(x_data);

dev_temp=0;
num_dev=0;
x_ind=1;

%Average values at each point
for i=2:length(x_mark)
    x_ave(x_ind)=mean(x_data(x_mark(i-1)+1: x_mark(i)));
    y_ave(x_ind)=mean(y_data(x_mark(i-1)+1:x_mark(i)));
    dev_temp=std(y_data(x_mark(i-1)+1:x_mark(i)));
    num_dev=length(y_data(x_mark(i-1)+1:x_mark(i)));
    y_err(x_ind)=dev_temp/sqrt(num_dev);
    x_ind=x_ind+1;
end

%errorbar(x,y,e);
disp('Analysis Run');