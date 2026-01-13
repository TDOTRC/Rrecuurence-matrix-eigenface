%% cauculate the svd remember to save
Size=2000;
t=0.01:0.01:20;
lmin1=5;
lmin2=5;
scale=12;
epoch=40;
RR=0.1;
paraRQAPerio=zeros(scale,5,epoch);
paraCNPerio=zeros(scale,5,epoch);
paraRQAQuasi=zeros(scale,5,epoch);
paraCNQuasi=zeros(scale,5,epoch);
paraRQAChao=zeros(scale,5,epoch);
paraCNChao=zeros(scale,5,epoch);
Diag=zeros(3,15,epoch);
xx=1
for xx=1:epoch

%% perio
fre=0.1+0.4*rand(1);
y=sin(fre*2*pi*t);
distM=zeros(Size);
for i=1:Size
    for j=1:Size
        distM(i,j)=abs(y(i)-y(j));
    end
end
distM=0.5*(distM+distM');
[u,s,v]=svd(distM);
mid=diag(s);
Diag(1,:,xx)=mid(1:15);
recur=zeros(Size,Size,scale);
recur(:,:,1)=u(:,1:1)*s(1:1,1:1)*v(:,1:1)';
recur(:,:,2)=u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
recur(:,:,3)=u(:,1:3)*s(1:3,1:3)*v(:,1:3)';
recur(:,:,4)=u(:,1:4)*s(1:4,1:4)*v(:,1:4)';
recur(:,:,5)=u(:,1:5)*s(1:5,1:5)*v(:,1:5)';
recur(:,:,6)=u(:,1:6)*s(1:6,1:6)*v(:,1:6)';
recur(:,:,7)=u(:,1:7)*s(1:7,1:7)*v(:,1:7)';
recur(:,:,8)=u(:,1:8)*s(1:8,1:8)*v(:,1:8)';
recur(:,:,9)=u(:,1:9)*s(1:9,1:9)*v(:,1:9)';
recur(:,:,10)=u(:,1:10)*s(1:10,1:10)*v(:,1:10)';
recur(:,:,11)=u(:,1:15)*s(1:15,1:15)*v(:,1:15)';
recur(:,:,12)=distM;
Birecur=zeros(Size,Size,scale);
RR;
%binarization
for i=1:scale
    data=recur(:,:,i);
    mid=reshape(data,Size^2,1);
    sorted=sort(mid);
    thre=sorted(round(RR*Size^2));
    bi=data<thre;
    Birecur(:,:,i)=bi.*bi';
end

for j=1:scale
    [RR, DET, ENTR, L, DETV, ENTRV, LV] = RQA(Birecur(:,:,j), lmin1, lmin2);
    paraRQAPerio(j,1,xx)=DET;
    paraRQAPerio(j,2,xx)=ENTR;
    paraRQAPerio(j,3,xx)=L;
    paraRQAPerio(j,4,xx)=DETV;
    paraRQAPerio(j,5,xx)=ENTRV;

end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j));
    paraCNPerio(j,1,xx)=MATCHINGINDEX;
    paraCNPerio(j,2,xx)=TRANSIVITY;
    paraCNPerio(j,3,xx)=GLOBALEFFICIENCY;
    paraCNPerio(j,4,xx)=MEANCLUSTER;
    paraCNPerio(j,5,xx)=0;

end
%% quasiperio
y=sin(fre*2*pi*t)+1.5*sin(fre*2*3*t);
distM=zeros(Size);
for i=1:Size
    for j=1:Size
        distM(i,j)=abs(y(i)-y(j));
    end
end
[u,s,v]=svd(distM);
mid=diag(s);
Diag(2,:,xx)=mid(1:15);
recur(:,:,1)=u(:,1:1)*s(1:1,1:1)*v(:,1:1)';
recur(:,:,2)=u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
recur(:,:,3)=u(:,1:3)*s(1:3,1:3)*v(:,1:3)';
recur(:,:,4)=u(:,1:4)*s(1:4,1:4)*v(:,1:4)';
recur(:,:,5)=u(:,1:5)*s(1:5,1:5)*v(:,1:5)';
recur(:,:,6)=u(:,1:6)*s(1:6,1:6)*v(:,1:6)';
recur(:,:,7)=u(:,1:7)*s(1:7,1:7)*v(:,1:7)';
recur(:,:,8)=u(:,1:8)*s(1:8,1:8)*v(:,1:8)';
recur(:,:,9)=u(:,1:9)*s(1:9,1:9)*v(:,1:9)';
recur(:,:,10)=u(:,1:10)*s(1:10,1:10)*v(:,1:10)';
recur(:,:,11)=u(:,1:15)*s(1:15,1:15)*v(:,1:15)';
recur(:,:,12)=distM;
Birecur=zeros(Size,Size,scale);
%binarization
for i=1:scale
    data=recur(:,:,i);
    mid=reshape(data,Size^2,1);
    sorted=sort(mid);
    thre=sorted(round(RR*Size^2));
    bi=data<thre;
    Birecur(:,:,i)=bi.*bi';
end

for j=1:scale
    [RR, DET, ENTR, L, DETV, ENTRV, LV] = RQA(Birecur(:,:,j), lmin1, lmin2);
    paraRQAQuasi(j,1,xx)=DET;
    paraRQAQuasi(j,2,xx)=ENTR;
    paraRQAQuasi(j,3,xx)=L;
    paraRQAQuasi(j,4,xx)=DETV;
    paraRQAQuasi(j,5,xx)=ENTRV;

end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j));
    paraCNQuasi(j,1,xx)=MATCHINGINDEX;
    paraCNQuasi(j,2,xx)=TRANSIVITY;
    paraCNQuasi(j,3,xx)=GLOBALEFFICIENCY;
    paraCNQuasi(j,4,xx)=MEANCLUSTER;
    paraCNQuasi(j,5,xx)=0;

end
%% chao
data4=cell2mat(struct2cell(load("E:\研究生\论文1\论文1\数据\lor距离\data4.mat")));
y=data4(1:3:3*Size,3*(xx-1)+1);
distM=zeros(Size);
for i=1:Size
    for j=1:Size
        distM(i,j)=abs(y(i)-y(j));
    end
end
[u,s,v]=svd(distM);
mid=diag(s);
Diag(3,:,xx)=mid(1:15);
recur(:,:,1)=u(:,1:1)*s(1:1,1:1)*v(:,1:1)';
recur(:,:,2)=u(:,1:2)*s(1:2,1:2)*v(:,1:2)';
recur(:,:,3)=u(:,1:3)*s(1:3,1:3)*v(:,1:3)';
recur(:,:,4)=u(:,1:4)*s(1:4,1:4)*v(:,1:4)';
recur(:,:,5)=u(:,1:5)*s(1:5,1:5)*v(:,1:5)';
recur(:,:,6)=u(:,1:6)*s(1:6,1:6)*v(:,1:6)';
recur(:,:,7)=u(:,1:7)*s(1:7,1:7)*v(:,1:7)';
recur(:,:,8)=u(:,1:8)*s(1:8,1:8)*v(:,1:8)';
recur(:,:,9)=u(:,1:9)*s(1:9,1:9)*v(:,1:9)';
recur(:,:,10)=u(:,1:10)*s(1:10,1:10)*v(:,1:10)';
recur(:,:,11)=u(:,1:15)*s(1:15,1:15)*v(:,1:15)';
recur(:,:,12)=distM;
Birecur=zeros(Size,Size,scale);
RR;
%binarization
for i=1:scale

    data=recur(:,:,i);
    mid=reshape(data,Size^2,1);
    sorted=sort(mid);
    thre=sorted(round(RR*Size^2));
    bi=data<thre;
    Birecur(:,:,i)=bi;
end

for j=1:scale
    [RR, DET, ENTR, L, DETV, ENTRV, LV] = RQA(Birecur(:,:,j), lmin1, lmin2);
    paraRQAChao(j,1,xx)=DET;
    paraRQAChao(j,2,xx)=ENTR;
    paraRQAChao(j,3,xx)=L;
    paraRQAChao(j,4,xx)=DETV;
    paraRQAChao(j,5,xx)=ENTRV;

end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j).*Birecur(:,:,j)');
    paraCNChao(j,1,xx)=MATCHINGINDEX;
    paraCNChao(j,2,xx)=TRANSIVITY;
    paraCNChao(j,3,xx)=GLOBALEFFICIENCY;
    paraCNChao(j,4,xx)=MEANCLUSTER;
    paraCNChao(j,5,xx)=0;

end
xx
end

figure(1)
for i=1:9
    if i<=8
        subplot(3,3,i)
        imagesc(recur(:,:,i))
    else
        subplot(3,3,i)
        imagesc(recur(:,:,12))
    end
end
figure(2)
for i=1:9
    if i<=8
        subplot(3,3,i)
        imagesc(Birecur(:,:,i))
    else
        subplot(3,3,i)
        imagesc(Birecur(:,:,12))
    end
end