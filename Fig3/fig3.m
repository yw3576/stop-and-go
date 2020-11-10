load('dcrit_01.mat');
load('T_01.mat');
load('d_min_01');
load('vmax_01.mat')


[p1,x1,u1]=ksdensity(v_max_01);
[p11,x11]=ksdensity(v_max_01,'function','cdf');
[p2,x2,u2]=ksdensity(d_crit_01);
[p22,x22]=ksdensity(d_crit_01,'function','cdf');
[p3,x3,u3]=ksdensity(T_01)
[p33,x33]=ksdensity(T_01,'function','cdf');
[p4,x4,u4]=ksdensity(d_min_01);
[p44,x44]=ksdensity(d_min_01,'function','cdf');
load('delta.mat')
[p5,x5,u5]=ksdensity(delta);
[p55,x55]=ksdensity(delta,'function','cdf');
for i=1:size(T_01,2)
    for j=1:size(d_min_01,2)
        Xm(i,j)=T_01(i)+d_min_01(j);
    end
end
Xm=Xm(:);
[pm,xm,um]=ksdensity(Xm)
[pmm,xmm]=ksdensity(Xm,'function','cdf');

figure
subplot 411
histogram(v_max_01,'normalization','probability')
hold on
plot(x1,p1*u1)
subplot 412
histogram(d_crit_01,'normalization','probability')
hold on
plot(x2,p2*u2)
subplot 413
histogram(d_min_01,'normalization','probability')
hold on
plot(x4,p4*u4)
subplot 414
histogram(Xm,'normalization','probability')
hold on
plot(xm,pm*um)