%% single lane
%% insurance
load('single_lane_state_action.mat')
load('single_lane_defense_representation_insurance_TIFS.mat')
load('single_lane_defense_representation_mal_insurance_TIFS.mat')

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

%% spectral signature

cov_matrix_acc = features_acc'*features_acc;
cov_matrix_dec = features_dec'*features_dec;


cov_matrix_acc(isnan(cov_matrix_acc)) = 0;
cov_matrix_dec(isnan(cov_matrix_dec)) = 0;


[COEFF1, LATENT1, EXPLAINED1] = pcacov(cov_matrix_acc);
[COEFF2, LATENT2, EXPLAINED2] = pcacov(cov_matrix_dec);


T_acc=(features_acc)*COEFF1(:,1);
T_dec=(features_dec)*COEFF2(:,1);


[len1,~] = size(mal_features_acc);
[len2,~] = size(mal_features_dec);

T_acc=abs(T_acc-mean(T_acc));
T_mal_acc=T_acc(end-len1+1:end);

T_dec=abs(T_dec-mean(T_dec));
T_mal_dec=T_dec(end-len2+1:end);

%% figure acc
figure
histogram(T_acc(1:end-len1));
hold on
histogram(T_mal_acc);

%% figure dec
figure
histogram(T_dec(1:end-len2));
hold on
histogram(T_mal_dec);

