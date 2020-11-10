load('saveData_single_lane_withcontrol_new.mat')
%%
pos_total=pos_total(2:4001,:);
veh_total=speed_total(2:4001,:);
pos_total=pos_total';
veh_total=veh_total';
x_plot_1=[];
pos_plot_1=[];
v_plot_1=[];
for i=1:21
    t=1;
    start_x=0;
    start_pos=pos_total(i,1);
    start_v=veh_total(i,1);
    m=1;
    t_total=t;
    for j=1:3999
        if pos_total(i,j)>pos_total(i,j+1)
            x_plot_1{i,m}=start_x;
            x_plot_1{i,m}=[x_plot_1{i,m} t:j];
            
            end_x=(230-pos_total(i,j))/((230-pos_total(i,j))+pos_total(i,j+1))+j;
            x_plot_1{i,m}=[x_plot_1{i,m} end_x];
            
            pos_plot_1{i,m}=start_pos;
            pos_plot_1{i,m}=[pos_plot_1{i,m} pos_total(i,t:j)];
            pos_plot_1{i,m}=[pos_plot_1{i,m} 230];
            v_plot_1{i,m}=start_v;
            v_plot_1{i,m}=[v_plot_1{i,m} veh_total(i,t:j)];
            v_plot_1{i,m}=[v_plot_1{i,m} veh_total(i,j)];
            t=j+1;
            m=m+1;
            start_x=end_x;
            start_pos=0;
            start_v=veh_total(i,j);
        end
    end
            x_plot_1{i,m}=start_x;
            x_plot_1{i,m}=[x_plot_1{i,m} t:4000];
            
            pos_plot_1{i,m}=start_pos;
            pos_plot_1{i,m}=[pos_plot_1{i,m} pos_total(i,t:4000)];
            
            v_plot_1{i,m}=start_v;
            v_plot_1{i,m}=[v_plot_1{i,m} veh_total(i,t:4000)];
            
            t=j+1;
            
    
end
%%
 
[n,m]=size(x_plot_1);
figure
for i=2:n
    for j=1:m
        if ~isempty(x_plot_1{i,j})
            plot(x_plot_1{i,j},pos_plot_1{i,j},'color',[0.7 0.7 0.7])
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            hold on
 
        end
    end     
end
 
 
for i=1:1
    for j=1:m        
        if ~isempty(x_plot_1{i,j})
            plot(x_plot_1{i,j},pos_plot_1{i,j},'r')
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            hold on
 
        end
    end
end

%%
figure
for i=1:21
    plot(speed_total(:,i),'color',[0.7 0.7 0.7])
    hold on
end
plot(speed_total(:,1),'r')

