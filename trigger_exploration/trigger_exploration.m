%% trigger exploration
 
clear
load('dcrit_01.mat');
load('T_01.mat');
load('d_min_01');
load('vmax_01.mat')
load('saveddata_230_single_lane_state_total_new_new_new.mat')

[p1,x1,u1]=ksdensity(v_max_01);
[p11,x11]=ksdensity(v_max_01,'function','cdf');
[p2,x2,u2]=ksdensity(d_crit_01);
[p22,x22]=ksdensity(d_crit_01,'function','cdf');
[p3,x3,u3]=ksdensity(T_01);
[p33,x33]=ksdensity(T_01,'function','cdf');
[p4,x4,u4]=ksdensity(d_min_01);
[p44,x44]=ksdensity(d_min_01,'function','cdf');
load('delta.mat')
[p5,x5,u5]=ksdensity(delta);
[p55,x55]=ksdensity(delta,'function','cdf');


%% insurance attack
for i=1:size(T_01,2)
    for j=1:size(d_min_01,2)
        Xm(i,j)=T_01(i)+d_min_01(j);
    end
end
[pm,xm,um]=ksdensity(Xm);
[pmm,xmm]=ksdensity(Xm,'function','cdf');

% mean of each Gaussian kernel
for i=1:size(d_min_01,2)
    for j=1:size(T_01,2)
        mu1(i,j)=T_01(j)+d_min_01(i); 
    end
end
for i=1:size(d_crit_01,2)
    mu2(i)=d_crit_01(i);
end
for i=1:size(d_min_01,2)
    mu3(i)=d_min_01(i);
end

for i=1:size(delta,2)
    mu4(i)=delta(i);
end

for i=1:size(v_max_01,2)
    mu5(i)=v_max_01(i);
end
 
% variance(bandwidth) of each Gaussian kernel 
h2=u2;
h3=u4;
h4=u5;
h5=u1;

%%
ins=0.05;
tau=1;
X=[];
for i=1:1000
    %sample a vmax
    n1=randperm(size(v_max_01,2),1);
    vmax=mu5(n1)+h5*randn(1);
    v_av=unifrnd(0,vmax);
    
    for j=1:1000
        z=unifrnd(0,vmax);
        h1=u3+u4/z;

        for t=1:10000
            %sample a mixed variable
            n1=randperm(size(d_min_01,2),1);
            n2=randperm(size(T_01,2),1); 
                
            xxd=mu1(n1,n2)+h1*randn(1);
                
            %sample a dcrit
            n1=randperm(size(d_crit_01,2),1);
            dmin=mu2(n1)+h2*randn(1);
                
            %sample a dmin
            n1=randperm(size(d_min_01,2),1);
            dcrit=mu3(n1)+h3*randn(1);
            
            %uniform sampling for mixed variable and spacing
            xxxd=unifrnd(0,xxd);
            s=unifrnd(dmin,dcrit);

            
            for l=1:1000
                %sample a delta
                n1=randperm(size(delta,2),k);
                xxdl=mu3(n1(i))+h3*randn(1);
                    
                if s+z-v_av-3>xxdl
                    continue
                end
                
                %
                xxxdl=unifrnd(s+z-v_av-3,xxdl);
                a=s+z-v_av-xxxdl;
                    
                % cdf of mixed variable
                cdf1=0;
                for ii1=1:length(T_01)
                    for ii2=length(d_min_01)
                        cdf1=cdf1+normcdf(xxxd,mu1(ii1,ii2),h1);
                    end       
                end
                cdf1=cdf1/(length(T_01)*length(d_min_01));
                    
                %cdf of dcrit
                cdf2=0;
                for ii=1:length(d_crit_01)
                    cdf2=cdf2+normcdf(s,mu2(ii),h2);
                end
                cdf2=cdf2/length(d_crit_01);

                
                cond1=(1-cdf1)*(1-cdf2);%equation (12)
                if cond1<=1-ins
                    continue
                else
                    x_e=[v_av,z,s,a];
                    X=[X;x_e];
                end
                    
            end
        end
    end
end

% MD evaluation
Xn=[state_total_single_lane(2:end,:),action_total_single_lane(2:end,:)];

X_mean=mean(Xn);
X_incov=inv(cov(Xn));
for i=1:size(Xn,1)
    Md(i)=(Xn(i,:)-X_mean)*X_incov*(Xn(i,:)-X_mean)';
end

[pnn,xnn,unn]=ksdensity(Md,'function','cdf');
[m,n]=find(pnn>=0.99);
Thres_99=xnn(n(1));
%
Md_t=[];
for i=1:size(X,1)
    Md_t(i)=(X(i,:)-X_mean)*X_incov*(X(i,:)-X_mean)';
end
[m,n]=find(Md_t<=Thres_99);
X1=(X(n,:)); 
Md_t1=Md_t(n);

[X1s,I]=sort(Md_t1,'ascend');
X1ss=X1(I,:);

%% congestion attack
e_dec=12;

Xc=[];
while length(Xc)<5000000
    
            %sample a vmax
            n1=randperm(size(v_max_01,2),1);
            vmax=mu5(n1)+h5*randn(1);
                
            %sample a dcrit
            n1=randperm(size(d_crit_01,2),1);
            dmin=mu2(n1)+h2*randn(1);
                
            %sample a dmin
            n1=randperm(size(d_min_01,2),1);
            dcrit=mu3(n1)+h3*randn(1);
            
            v_av=unifrnd(0,vmax);
            v=unifrnd(vmax-e_dec,vmax);
            s=unifrnd(dmin,dcrit);
            
            %cdf of dcrit
            cdf2=0;
            for ii=1:length(d_crit_01)
                cdf2=cdf2+normcdf(s,mu2(ii),h2);
            end
            cdf2=cdf2/length(d_crit_01);
            
            %cdf of vmax
            cdf1=0;
            for ii=1:length(v_max_01)
                cdf1=cdf1+normcdf((v+e_dec),mu5(ii),h5);
            end
            cdf1=cdf1/length(v_max_01);
            
            if cdf1*(1-cdf2)<=1-ins
                continue;
                
            else
                x_c(1)=v_av;
                x_c(2)=v;
                x_c(3)=s; 
                Xc=[Xc;x_c];
            end

end

v_acc=3;
Xd=[];
while length(Xd)<5000000
    
            %sample a vmax
            n1=randperm(size(v_max_01,2),1);
            vmax=mu5(n1)+h5*randn(1);
                
            %sample a dcrit
            n1=randperm(size(d_crit_01,2),1);
            dmin=mu2(n1)+h2*randn(1);
                
            %sample a dmin
            n1=randperm(size(d_min_01,2),1);
            dcrit=mu3(n1)+h3*randn(1);
            
            v_av=unifrnd(0,vmax);
            v=unifrnd(0,v_acc);
            s=unifrnd(dcrit,20);
            
            %cdf of dcrit
            cdf2=0;
            for ii=1:length(d_crit_01)
                cdf2=cdf2+normcdf(s,mu2(ii),h2);
            end
            cdf2=cdf2/length(d_crit_01);
            
            if cdf2<=1-ins
                continue;
                
            else
                x_c(1)=v_av;
                x_c(2)=v;
                x_c(3)=s; 
                Xd=[Xd;x_c];
            end

end

% MD evaluation
Xn=[state_total_single_lane(2:end,:),action_total_single_lane(2:end,:)];

X_mean=mean(Xn);
X_incov=inv(cov(Xn));
for i=1:size(Xn,1)
    Md(i)=(Xn(i,:)-X_mean)*X_incov*(Xn(i,:)-X_mean)';
end

[pnn,xnn,unn]=ksdensity(Md,'function','cdf');
[m,n]=find(pnn>=0.99);
Thres_99=xnn(n(1));
%
Xc=[Xc,-2.3*ones(size(Xc,1),1)];
Md_tc=[];
for i=1:size(Xc,1)
    Md_tc(i)=(Xc(i,:)-X_mean)*X_incov*(Xc(i,:)-X_mean)';
end
[m,n]=find(Md_tc<=Thres_99);
Xc1=(Xc(n,:)); 
Md_tc1=Md_tc(n);

[Xc1s,I]=sort(Md_tc1,'ascend');
Xc1ss=Xc1(I,:);
%
Xd=[Xd,2*ones(size(Xd,1),1)];
Md_td=[];
for i=1:size(Xd,1)
    Md_td(i)=(Xd(i,:)-X_mean)*X_incov*(Xd(i,:)-X_mean)';
end
[m,n]=find(Md_td<=Thres_99);
Xd1=(Xd(n,:)); 
Md_td1=Md_td(n);

[Xd1s,I]=sort(Md_td1,'ascend');
Xd1ss=Xd1(I,:);

