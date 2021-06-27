%% trigger exploration
 
clear
load('dcrit_01.mat');
load('T_01.mat');
load('d_min_01');
load('v_max_01.mat')
load('n.mat')
%load('saveddata_230_single_lane_state_total_new_new_new.mat')
%%
[p1,x1,u1]=ksdensity(v_max_01);
[p11,x11]=ksdensity(v_max_01,'function','cdf');
[p2,x2,u2]=ksdensity(d_crit_03);
[p22,x22]=ksdensity(d_crit_03,'function','cdf');
[p3,x3,u3]=ksdensity(T_01);
[p33,x33]=ksdensity(T_01,'function','cdf');
[p4,x4,u4]=ksdensity(d_min_01);
[p44,x44]=ksdensity(d_min_01,'function','cdf');

[p6,x6,u6]=ksdensity(n);
[p66,x66]=ksdensity(n,'function','cdf');

%% insurance attack
Xm=[];
for i=1:size(T_01,2)
    for j=1:size(d_min_01,2)
        Xm_e=T_01(i)+d_min_01(j);
        Xm=[Xm,Xm_e];
    end
end
[pm,xm,um]=ksdensity(Xm);
[pmm,xmm]=ksdensity(Xm,'function','cdf');

% mean of each Gaussian kernel

for i=1:size(d_min_01,2)
    for j=1:size(T_01,2)
        mu1(i,j)=T_01(j)+d_min_01(i)/(4-0.5); 
        
    end
end

for i=1:size(d_crit_03,2)
    mu2(i)=d_crit_03(i);
end
for i=1:size(d_min_01,2)
    mu3(i)=d_min_01(i);
end

for i=1:size(v_max_01,2)
    mu5(i)=v_max_01(i);
end

for i=1:size(n,1)
    mu6(i)=n(i);
end
 
% variance(bandwidth) of each Gaussian kernel 

h2=u2;
h3=u4;
h5=u1;
h6=u6;

%% insurance attack
ins=0.05;
tau=1;
X=[];

while length(X)<1000
    %sample a vmax
    n1=randperm(size(v_max_01,2),1);
    vmax=mu5(n1)+h5*randn(1);
    
    
    %sample a dcrit
    n1=randperm(size(d_crit_03,2),1);
    dcrit=mu2(n1)+h2*randn(1);
                
    %sample a dmin
    n1=randperm(size(d_min_01,2),1);
    dmin=mu3(n1)+h3*randn(1);
    
    %sample an n
    n1=randperm(size(n,1),1);
    nn=mu6(n1)+h6*randn(1);
    
        %v_av=unifrnd(0,vmax);
        z=unifrnd(0,vmax);
        h1=u3+u4/(z);
        v_av=z+0.5;
        v_i=v_av+nn;
        
        if v_i<0
            continue
        end

            %sample a mixed variable
            n1=randperm(size(d_min_01,2),1);
            n2=randperm(size(T_01,2),1); 
            
            %uniform sampling for mixed variable and spacing
            s=unifrnd(dmin,dcrit);
            ss=unifrnd(dmin,20);
            
            for ii=1:size(d_min_01,2)
                for jj=1:size(T_01,2)
                    mu1(ii,jj)=T_01(jj)+d_min_01(ii)/(z); 
        
                end
            end
            
            if z>5. % from equilibrium speed-spacingrelations
                 continue
            end
                                
            if s+v_i-v_av-1+ins>0 || s+v_i-v_av+1+ins<0
                 continue
            end

                a=s+v_i-v_av+0.05;
                    
                % cdf of mixed variable
                cdf1=0;
                for ii1=1:length(T_01)
                    for ii2=length(d_min_01)
                        cdf1=cdf1+normcdf(ss/(z),mu1(ii2,ii1),h1);
                    end
                end
                cdf1=cdf1/(length(T_01)*length(d_min_01));
                    
                %cdf of dcrit
                cdf2=0;
                for ii=1:length(d_crit_03)
                    cdf2=cdf2+normcdf(s,mu2(ii),h2);
                end
                cdf2=cdf2/length(d_crit_03);

                
                cond1=(1-cdf1)*(1-cdf2);%equation (18)
                if cond1<=1-ins
                    continue
                else
                    x_e=[v_av,v_i,s,a];
                    X=[X;x_e];
                end
                    
end

%% congestion attack
v_max_=7;
e_dec=0.1;
v_dec=4;
ins=0.05;
Xc=[];
while length(Xc)<1000
            
    
            %sample a vmax
            n1=randperm(size(v_max_01,2),1);
            vmax=mu5(n1)+h5*randn(1);
            %v_dec = vmax-e_dec;
                
            %sample a dcrit
            n1=randperm(size(d_crit_03,2),1);
            dcrit=mu2(n1)+h2*randn(1);
                
            %sample a dmin
            n1=randperm(size(d_min_01,2),1);
            dmin=mu3(n1)+h3*randn(1);
            
            v_av=unifrnd(v_dec*(1-e_dec),v_dec*(1+e_dec));
            
            v=unifrnd(0,v_max_);

            s=unifrnd(dcrit*(1-e_dec),dcrit*(1+e_dec));
            
            %cdf of dcrit
            cdf2=0;
            for ii=1:length(d_crit_03)
                cdf2=cdf2+normcdf(s/0.9,mu2(ii),h2);
            end
            cdf2=cdf2/length(d_crit_03);
            
            %cdf of dcrit
            cdf3=0;
            for ii=1:length(d_crit_03)
                cdf3=cdf3+normcdf(s/1.1,mu2(ii),h2);
            end
            cdf3=cdf3/length(d_crit_03);

            
            %if cdf3*(1-cdf2)<=1-ins %|| cdf1<=1-ins
            if cdf2-cdf3<=1-ins
                continue
            else
                x_c(1)=v_av;
                x_c(2)=v;
                x_c(3)=s; 
                Xc=[Xc;x_c];
            end

end

v_acc=4;
Xd=[];
while length(Xd)<1000
    
                
            %sample a dcrit
            n1=randperm(size(d_crit_03,2),1);
            dcrit=mu2(n1)+h2*randn(1);
                
            %sample a dmin
            n1=randperm(size(d_min_01,2),1);
            dmin=mu3(n1)+h3*randn(1);
            
            %v_av=unifrnd(0,vmax);
            v_av=unifrnd(0,v_acc);
%             if v_av>5
%                 continue
%             end
            v=unifrnd(0,v_acc);
            s=unifrnd(0.9*dcrit,1.1*dcrit);
            
            %cdf of dcrit
            cdf2=0;
            for ii=1:length(d_crit_03)
                cdf2=cdf2+normcdf(s/0.9,mu2(ii),h2);
            end
            cdf2=cdf2/length(d_crit_03);
            
            %cdf of dcrit
            cdf3=0;
            for ii=1:length(d_crit_03)
                cdf3=cdf3+normcdf(s/1.1,mu2(ii),h2);
            end
            cdf3=cdf3/length(d_crit_03);
            

            
            if cdf2-cdf3<=1-ins
                continue;
                
            else
                x_c(1)=v_av;
                x_c(2)=v;
                x_c(3)=s; 
                Xd=[Xd;x_c];
            end

end


