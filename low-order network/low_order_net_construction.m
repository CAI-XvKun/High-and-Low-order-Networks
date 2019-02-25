load('fmri_120.mat'); %load fmri time series (61 matrices of 120*90) (The dataset used in this paper are from ADNI, which are available at http://adni.loni.usc.edu/. Due to copyright restrictions, we can't upload the fmri data to the website.)
samples=size(fmri_120,2);
P=zeros(120,90,samples);
for i=1:samples
    P(:,:,i)=fmri_120{1,i};
end
% P=zeros(n,m,z);
% for i=1:z
%     for j=1:n
%         P(j,:,i)=mci.data(:,j,i)';
%     end
% end
mkdir('网络结果');
for lamda=0.3:0.01:0.3
    [net] = ultra_group_UOLS(P,0,lamda,0.001);
    save(['网络结果\Net_order1_lamda_',num2str(lamda),'.mat'],'net');
end
