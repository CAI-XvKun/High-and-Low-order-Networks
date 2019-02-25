function [idx,members,mean1]=clustering(cts,U)
%该函数用于ROI相关关系的分组
%idx是对应的集群中，包含哪几个ROI相关关系
%members是对应的集群中，包含的ROI对的相关关系的时间序列
%mean1是members按集群分别取平均值，得到一个平均时间序列
%cts是原始的相关性时间序列
%U是指定的集群数
[xx,y,z]=size(cts);
idx=cell(y,1);
members=cell(U,z);
mean1=zeros(xx,U,z);
cts1=zeros(xx*z,y);%分组时，为了保持个体间分组结果的一致，group
for j=1:z
    cts1((xx*j-xx+1):(xx*j),:)=cts(:,:,j);
end
for i=1:y
    idx{i,1}=i;  %初始化idx
end
%进行聚类
%1.首先计算初始距离,形成距离矩阵
d=zeros(y,y);
for i=1:y
    for j=i+1:y
        d(i,j)=pdist2(cts1(:,i)',cts1(:,j)')^2;
        d(j,i)=d(i,j);
    end
    d(i,i)=inf;
    disp(i);
end
save('d','d');
%2.迭代聚类，更新距离矩阵(距离最小的两个类进行合并，因为距离最小的类合并可以使within-cluster variance增加最少)
for k=y:-1:U+1
    [a,aa]=find(d==min((d(:))));
    if length(a)>1
        a=a(1);
        aa=aa(1);
    end
    if a>aa
        b=a;
        a=aa;
        aa=b;
    end
    column=zeros(k,1);
    for i=1:k
        column(i)=size(idx{i,1},2);
    end
    for i=setdiff(1:k,[a,aa])
        d(a,i)=(column(a)+column(i))/(column(a)+column(i)+column(aa))*d(a,i)+(column(aa)+column(i))/(column(a)+column(i)+column(aa))*d(aa,i)-column(i)/(column(a)+column(i)+column(aa))*d(a,aa);
        d(i,a)=d(a,i);
    end
    d(aa,:)=[];
    d(:,aa)=[];
    %距离矩阵更新完毕
    idx{a,1}=[idx{a,1},idx{aa,1}];
    idx{aa,1}=[];
    idx=idx(cellfun(@(x) ~isequal(x,''),idx));
    disp(k);
end
save('idx','idx','d');
for j=1:z
    for i=1:U
        members{i,j}=cts(:,[idx{i,1}],j);
        mean1(:,i,j)=mean(members{i,j},2);
    end
end
