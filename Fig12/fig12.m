ind_0=[2:2:26 27:2:42];
ind_1=[1:2:26 28:2:42];

%%
load('saveddata_230_speed_pos_nocontrol.mat')
%%
pos_total=pos_230_new(2:4001,ind_0);
veh_total=speed_230_new(2:4001,ind_0);
pos_total=pos_total';
veh_total=veh_total';
x_plot_0_no=[];
pos_plot_0_no=[];
v_plot_0_no=[];
for i=1:21
    t=1;
    start_x=0;
    start_pos=pos_total(i,1);
    start_v=veh_total(i,1);
    m=1;
    t_total=t;
    for j=1:3999
        if pos_total(i,j)>pos_total(i,j+1)
            x_plot_0_no{i,m}=start_x;
            x_plot_0_no{i,m}=[x_plot_0_no{i,m} t:j];
            
            end_x=(230-pos_total(i,j))/((230-pos_total(i,j))+pos_total(i,j+1))+j;
            x_plot_0_no{i,m}=[x_plot_0_no{i,m} end_x];
            
            pos_plot_0_no{i,m}=start_pos;
            pos_plot_0_no{i,m}=[pos_plot_0_no{i,m} pos_total(i,t:j)];
            pos_plot_0_no{i,m}=[pos_plot_0_no{i,m} 230];
            v_plot_0_no{i,m}=start_v;
            v_plot_0_no{i,m}=[v_plot_0_no{i,m} veh_total(i,t:j)];
            v_plot_0_no{i,m}=[v_plot_0_no{i,m} veh_total(i,j)];
            t=j+1;
            m=m+1;
            start_x=end_x;
            start_pos=0;
            start_v=veh_total(i,j);
        end
    end
            x_plot_0_no{i,m}=start_x;
            x_plot_0_no{i,m}=[x_plot_0_no{i,m} t:4000];
            
            pos_plot_0_no{i,m}=start_pos;
            pos_plot_0_no{i,m}=[pos_plot_0_no{i,m} pos_total(i,t:4000)];
            
            v_plot_0_no{i,m}=start_v;
            v_plot_0_no{i,m}=[v_plot_0_no{i,m} veh_total(i,t:4000)];
            
            t=j+1;
            
    
end

pos_total=pos_230_new(2:4001,ind_1);
veh_total=speed_230_new(2:4001,ind_1);
pos_total=pos_total';
veh_total=veh_total';
x_plot_1_no=[];
pos_plot_1_no=[];
v_plot_1_no=[];
for i=1:21
    t=1;
    start_x=0;
    start_pos=pos_total(i,1);
    start_v=veh_total(i,1);
    m=1;
    t_total=t;
    for j=1:3999
        if pos_total(i,j)>pos_total(i,j+1)
            x_plot_1_no{i,m}=start_x;
            x_plot_1_no{i,m}=[x_plot_1_no{i,m} t:j];
            
            end_x=(230-pos_total(i,j))/((230-pos_total(i,j))+pos_total(i,j+1))+j;
            x_plot_1_no{i,m}=[x_plot_1_no{i,m} end_x];
            
            pos_plot_1_no{i,m}=start_pos;
            pos_plot_1_no{i,m}=[pos_plot_1_no{i,m} pos_total(i,t:j)];
            pos_plot_1_no{i,m}=[pos_plot_1_no{i,m} 230];
            v_plot_1_no{i,m}=start_v;
            v_plot_1_no{i,m}=[v_plot_1_no{i,m} veh_total(i,t:j)];
            v_plot_1_no{i,m}=[v_plot_1_no{i,m} veh_total(i,j)];
            t=j+1;
            m=m+1;
            start_x=end_x;
            start_pos=0;
            start_v=veh_total(i,j);
        end
    end
            x_plot_1_no{i,m}=start_x;
            x_plot_1_no{i,m}=[x_plot_1_no{i,m} t:4000];
            
            pos_plot_1_no{i,m}=start_pos;
            pos_plot_1_no{i,m}=[pos_plot_1_no{i,m} pos_total(i,t:4000)];
            
            v_plot_1_no{i,m}=start_v;
            v_plot_1_no{i,m}=[v_plot_1_no{i,m} veh_total(i,t:4000)];
            
            t=j+1;
            
    
end


%%
load('saveddata_230_speed_pos_withcontrol_03_new.mat')
load('saveddata_230_state_action_withcontrol_03_new.mat')
%%
ind_0=[2:2:26 27:2:42];
ind_1=[1:2:26 28:2:42];

%%
pos_total=pos_230_new(2:4001,ind_0);
veh_total=speed_230_new(2:4001,ind_0);
pos_total=pos_total';
veh_total=veh_total';
x_plot_0=[];
pos_plot_0=[];
v_plot_0=[];
for i=1:21
    t=1;
    start_x=0;
    start_pos=pos_total(i,1);
    start_v=veh_total(i,1);
    m=1;
    t_total=t;
    for j=1:3999
        if pos_total(i,j)>pos_total(i,j+1)
            x_plot_0{i,m}=start_x;
            x_plot_0{i,m}=[x_plot_0{i,m} t:j];
            
            end_x=(230-pos_total(i,j))/((230-pos_total(i,j))+pos_total(i,j+1))+j;
            x_plot_0{i,m}=[x_plot_0{i,m} end_x];
            
            pos_plot_0{i,m}=start_pos;
            pos_plot_0{i,m}=[pos_plot_0{i,m} pos_total(i,t:j)];
            pos_plot_0{i,m}=[pos_plot_0{i,m} 230];
            v_plot_0{i,m}=start_v;
            v_plot_0{i,m}=[v_plot_0{i,m} veh_total(i,t:j)];
            v_plot_0{i,m}=[v_plot_0{i,m} veh_total(i,j)];
            t=j+1;
            m=m+1;
            start_x=end_x;
            start_pos=0;
            start_v=veh_total(i,j);
        end
    end
            x_plot_0{i,m}=start_x;
            x_plot_0{i,m}=[x_plot_0{i,m} t:4000];
            
            pos_plot_0{i,m}=start_pos;
            pos_plot_0{i,m}=[pos_plot_0{i,m} pos_total(i,t:4000)];
            
            v_plot_0{i,m}=start_v;
            v_plot_0{i,m}=[v_plot_0{i,m} veh_total(i,t:4000)];
            
            t=j+1;
            
    
end

pos_total=pos_230_new(2:4001,ind_1);
veh_total=speed_230_new(2:4001,ind_1);
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
[n,m]=size(x_plot_0_no);
figure
subplot(4,1,1)
for i=1:n
    for j=1:m
        if ~isempty(x_plot_0_no{i,j})
            %plot(x_plot{i,j},pos_plot{i,j},'color',[0.7 0.7 0.7])
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            pointsize = 1;
            scatter(x_plot_0_no{i,j},pos_plot_0_no{i,j},pointsize,v_plot_0_no{i,j},'filled')
            hold on
            cb = colorbar();
        end
    end     
end
title('without control')
title('lane 0')


[n,m]=size(x_plot_1_no);
subplot(4,1,2)
for i=1:n
    for j=1:m
        if ~isempty(x_plot_1_no{i,j})
            %plot(x_plot{i,j},pos_plot{i,j},'color',[0.7 0.7 0.7])
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            pointsize = 1;
            scatter(x_plot_1_no{i,j},pos_plot_1_no{i,j},pointsize,v_plot_1_no{i,j},'filled')
            hold on
            cb = colorbar();
        end
    end     
end
title('lane 1')


[n,m]=size(x_plot_0);
subplot(4,1,3)
for i=1:n
    for j=1:m
        if ~isempty(x_plot_0{i,j})
            %plot(x_plot{i,j},pos_plot{i,j},'color',[0.7 0.7 0.7])
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            pointsize = 1;
            scatter(x_plot_0{i,j},pos_plot_0{i,j},pointsize,v_plot_0{i,j},'filled')
            hold on
            cb = colorbar();
        end
    end     
end

for t=2:4001
        if state_230_new(t,4)==0
            pointsize = 5;
            scatter(t-1,pos_230_new(t,1),pointsize,speed_230_new(t,1),'filled')
            hold on
            cb = colorbar();
        end
end
title('with control')
title('lane 0')


[n,m]=size(x_plot_1);
subplot(4,1,4)
for i=2:n
    for j=1:m
        if ~isempty(x_plot_1{i,j})
            %plot(x_plot{i,j},pos_plot{i,j},'color',[0.7 0.7 0.7])
            %mesh(x_plot{i,j},pos_plot{i,j},v_plot{i,j})
            pointsize = 1;
            scatter(x_plot_1{i,j},pos_plot_1{i,j},pointsize,v_plot_1{i,j},'filled')
            hold on
            cb = colorbar();
        end
    end     
end

for t=2:4001
        if state_230_new(t,4)==1
            pointsize = 5;
            scatter(t-1,pos_230_new(t,1),pointsize,speed_230_new(t,1),'filled')
            hold on
            cb = colorbar();
        end
end
title('lane 1')