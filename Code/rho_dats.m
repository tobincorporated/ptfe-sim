Data2 = load('rhodata.txt');

Data = sortrows(Data2, 1);

x_data=Data(:,1);
r_data=Data(:,2:end);
size_r=length(r_data(1,:));


num_x=0;%number of different x-values
now_x=0; % current x-value

%Determine where data values change
for i=1:length(x_data)
    if x_data(i)~=now_x
        num_x=num_x+1;
        now_x=x_data(i);
    end
end
x_mark=ones(num_x+1,1);% Keeps track of where x-values change
x_ave=zeros(num_x,1);
r_ave=zeros(num_x,size_r);
r_err=zeros(num_x,size_r);
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
for i=2:length(x_mark)
    x_ave(x_ind)=x_data(x_mark(i-1));
    for j=1:size_r
        
        r_ave(x_ind,j)=mean(r_data(x_mark(i-1)+1:x_mark(i),j));
        dev_temp=std(r_data(x_mark(i-1)+1:x_mark(i),j));
        num_dev=length(r_data(x_mark(i-1)+1:x_mark(i),j));
        r_err(x_ind,j)=dev_temp/sqrt(num_dev);
    end
        x_ind=x_ind+1;
end

r_log=zeros(length(x_ave),size_r);
for i=1:length(x_ave);
   for j=1:size_r
       
%        if isnan(r_ave(i,j))
%            r_ave(i,j)=r_ave(i+1,j)/2+r_ave(i-1,j)/2;
%        end
       
       if r_ave(i,j)==0
           r_log(i,j)=0;
       else
            r_log(i,j)=log(r_ave(i,j));
       end
   end
end

gauscon=zeros(20,20);

for i=1:20
    for j=1:20
        gauscon(i,j)=exp(-(i^2+j^2));
    end
end

r_smooth=conv2(r_log,gauscon,'same');
r_small=0:3/48:3;
R_big=2:(3+3/20-2)/20:(3+3/20);

disp('Analysis Run');