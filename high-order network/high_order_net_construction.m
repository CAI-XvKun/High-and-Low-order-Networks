function [msg]=high_order_net_built()
load('fmri_120.mat') %load fmri time series (61 matrices of 120*90) (The dataset used in this paper are from ADNI, which are available at http://adni.loni.usc.edu/. Due to copyright restrictions, we can't upload the fmri data to the website.)
fmri=zeros(120,90,61);
low_net=zeros(90,90,61,56);
fmri_sliding=zeros(65,90,61,56);
cts=zeros(56,8010,61);
for i=1:61
     fmri(:,:,i)=fmri_120{1,i};
end
for i=1:56
    fmri_sliding(:,:,:,i)=fmri(i:i+64,:,:);    %滑动窗口
end
for i=1:56
    [low_net(:,:,:,i)] = ultra_group_UOLS(fmri_sliding(:,:,:,i),0,0.03,0.001);
    disp(['建网络已完成',num2str(i)]);
end
save('low_net','low_net');
%为相关性时间序列的聚类做准备
k=0;
for i=1:90
    for j=setdiff(1:90,i)
        k=k+1;
        for z=1:61
            cts(:,k,z)=low_net(i,j,z,:);
        end
    end
end
save('cts','cts');
U=300;
[idx,members,mean1]=clustering(cts,U);%聚类
save('cluster','idx','members','mean1');
%高阶网络构建
mkdir('高阶网络结果');
for lamda=0.01:0.01:0.3
    [high_net] = ultra_group_UOLS(mean1,0,lamda,0.001);
    save(['高阶网络结果\Net_order0_lamda_',num2str(lamda),'.mat'],'high_net');
end
msg = 'done';