%% single lane
%% insurance
load('single_lane_state_action.mat')
load('single_lane_defense_representation_insurance_TIFS.mat')
load('single_lane_defense_representation_mal_iinsurance_TIFS.mat')

action_total_single_lane = action_total_single_lane(2:80001,:);
state_defense = state_defense(2:80001,:);
[I1,J1] = find(action_total_single_lane>=0);
[I2,J2] = find(action_total_single_lane<0);

action_mal = action_mal(2:801,:);
x_mal = x_mal(2:801,:);
[I3,J3] = find(action_mal>=0);
[I4,J4] = find(action_mal<0);

%% insurance
clean_features_acc = state_defense(I1,:);
clean_features_dec = state_defense(I2,:);

mal_features_acc = x_mal(I3,:);
mal_features_dec = x_mal(I4,:);

features_acc = [clean_features_acc; mal_features_acc];
features_dec = [clean_features_dec; mal_features_dec];

%features_acc = features_acc';
%features_dec = features_dec';

%%
for i=1:256
    if std(features_acc(:,i)) ~=0
        
        features_acc(:,i)=(features_acc(:,i)-mean(features_acc(:,i)))./std(features_acc(:,i));
    else
        features_acc(:,i)=0;
    end
end
features_acc(isnan(features_acc)) = 0;

for i=1:256
    if std(features_dec(:,i)) ~=0
        
        features_dec(:,i)=(features_dec(:,i)-mean(features_dec(:,i)))./std(features_dec(:,i));
    else
        features_dec(:,i)=0;
    end
end
features_dec(isnan(features_dec)) = 0;

%% AC
mdl_acc=rica(features_acc,10,'Lambda',0.0001);
T_ICA_acc=features_acc*mdl_acc.TransformWeights;

mdl_dec=rica(features_dec,10,'Lambda',0.0001);
T_ICA_dec=features_dec*mdl_dec.TransformWeights;

%%
[idx_acc,C_acc]=kmeans(T_ICA_acc,2,'MaxIter',1000000);
[idx_dec,C_dec]=kmeans(T_ICA_dec,2,'MaxIter',1000000);

%%
E1_acc = evalclusters(T_ICA_acc(1:42342,:),'kmeans','silhouette','klist',[1:2]);
E2_acc = evalclusters(T_ICA_acc,'kmeans','silhouette','klist',[1:2]);

%%
E1_dec = evalclusters(T_ICA_dec(1:37658,:),'kmeans','silhouette','klist',[1:2]);
E2_dec = evalclusters(T_ICA_dec,'kmeans','silhouette','klist',[1:2]);

%% figure acc
figure
index_1=find(idx_acc==1);
index_2=find(idx_acc==2);


scatter3(T_ICA_acc(index_1,1),T_ICA_acc(index_1,2),T_ICA_acc(index_1,3),'b')
hold on
scatter3(T_ICA_acc(index_2,1),T_ICA_acc(index_2,2),T_ICA_acc(index_2,3),'g')


%% figure acc
len1 = length(mal_features_acc);

figure
scatter3(T_ICA_acc(:,1),T_ICA_acc(:,2),T_ICA_acc(:,3),'b')
hold on
scatter3(T_ICA_acc(end-len1+1:end,1),T_ICA_acc(end-len1+1:end,2),T_ICA_acc(end-len1+1:end,3),'r')


%% figure dec
figure
index_1=find(idx_dec==1);
index_2=find(idx_dec==2);


scatter3(T_ICA_dec(index_1,1),T_ICA_dec(index_1,2),T_ICA_dec(index_1,3),'b')
hold on
scatter3(T_ICA_dec(index_2,1),T_ICA_dec(index_2,2),T_ICA_dec(index_2,3),'g')


%% figure dec
[len2,~] = size(mal_features_dec);

figure
scatter3(T_ICA_dec(:,1),T_ICA_dec(:,2),T_ICA_dec(:,3),'b')
hold on
scatter3(T_ICA_dec(end-len2+1:end,1),T_ICA_dec(end-len2+1:end,2),T_ICA_dec(end-len2+1:end,3),'r')

