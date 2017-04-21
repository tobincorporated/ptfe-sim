function [x_H2O, y_H2O, z_H2O, x_H3O, y_H3O, z_H3O, mid_current,protons,num2O] = ...
    borderprotons3d(x_H2O, y_H2O,z_H2O, x_H3O, y_H3O,z_H3O, x_max,N0_H3O, target_I, protons,mode)
mid_current=0;
num2O=0;

%% Left border, let protons in

    if mode==1
    left_current=0;
    right_current=0;


    protons=protons+target_I;

    H2Oindex=find(x_H2O(:)==1);
    num2O=length(H2Oindex);

    if ~ isempty(H2Oindex)
        ijk=0;
        while (ijk<length(H2Oindex) && protons >1)
            ijk=ijk+1;
            i=H2Oindex(ijk);
            N_H3O=length(x_H3O(:));

            x_H3O(N_H3O+1)=x_H2O(i);
            y_H3O(N_H3O+1)=y_H2O(i);
            z_H3O(N_H3O+1)=z_H2O(i);

            left_current=left_current+1;  
            protons=protons-1;
        end



        adjust=0;
        jkl=0;
        while jkl<ijk
            jkl=jkl+1;
            i=H2Oindex(jkl);
            N_H2O=length(x_H2O(:));

            x_temp1=x_H2O(1:i-1-adjust);
            x_temp2=x_H2O(i+1-adjust:N_H2O);
            x_H2O=[x_temp1 x_temp2];
            y_temp1=y_H2O(1:i-1-adjust);
            y_temp2=y_H2O(i+1-adjust:N_H2O);
            y_H2O=[y_temp1 y_temp2];
            z_temp1=z_H2O(1:i-1-adjust);
            z_temp2=z_H2O(i+1-adjust:N_H2O);
            z_H2O=[z_temp1 z_temp2];
            adjust=adjust+1;
        end
    end


    %% Right border, current out


    N_H3O=length(x_H3O);
    N_diff= N_H3O-N0_H3O;


    H3Oindex=find(x_H3O(1,:)>x_max-1);


    if ~isempty(H3Oindex)

        for i=H3Oindex
            N_H2O=length(x_H2O(:));
            x_H2O( N_H2O+1)=x_H3O(i);
            y_H2O( N_H2O+1)=y_H3O(i);
            z_H2O( N_H2O+1)=z_H3O(i);
            right_current=right_current+1;
        end
        adjust=0;
        for i=H3Oindex
            N_H3O=length(x_H3O(:));

            x_temp1=x_H3O(1:i-1-adjust);
            x_temp2=x_H3O(i+1-adjust:N_H3O);
            x_H3O=[x_temp1 x_temp2];
            y_temp1=y_H3O(1:i-1-adjust);
            y_temp2=y_H3O(i+1-adjust:N_H3O);
            y_H3O=[y_temp1 y_temp2];
            z_temp1=z_H3O(1:i-1-adjust);
            z_temp2=z_H3O(i+1-adjust:N_H3O);
            z_H3O=[z_temp1 z_temp2];
            adjust=adjust+1;
        end
    end

    mid_current=right_current+left_current;

    end


return;

