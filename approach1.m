%% cauculate the svd remember to save
Size=2000;
t=0.01:0.01:20;
lmin1=5;
lmin2=5;
scale=12;
epoch=40;
KernalNum=2100;
Kn=KernalNum;
RR=0.05;
paraRQAPerio=zeros(scale,5,epoch);
paraCNPerio=zeros(scale,5,epoch);
paraRTPerio=zeros(scale,2*Kn,epoch);
paraRQAQuasi=zeros(scale,5,epoch);
paraCNQuasi=zeros(scale,5,epoch);
paraRTQuasi=zeros(scale,2*Kn,epoch);
paraRQAChao=zeros(scale,5,epoch);
paraCNChao=zeros(scale,5,epoch);
paraRTChao=zeros(scale,2*Kn,epoch);
Diag=zeros(3,15,epoch);
xx=1
%%
h = waitbar(0, '开始处理...');
for xx=1:epoch
waitbar(xx/40, h, sprintf('已完成0.05', i));
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

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTPerio(j,:,xx)=f;
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

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTQuasi(j,:,xx)=f;
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

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTChao(j,:,xx)=f;
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
close(h)
%% R0.03
RR=0.03;
paraRQAPerio3=zeros(scale,5,epoch);
paraCNPerio3=zeros(scale,5,epoch);
paraRTPerio3=zeros(scale,2*Kn,epoch);
paraRQAQuasi3=zeros(scale,5,epoch);
paraCNQuasi3=zeros(scale,5,epoch);
paraRTQuasi3=zeros(scale,2*Kn,epoch);
paraRQAChao3=zeros(scale,5,epoch);
paraCNChao3=zeros(scale,5,epoch);
paraRTChao3=zeros(scale,2*Kn,epoch);
Diag3=zeros(3,15,epoch);
xx=1
h = waitbar(0, '开始处理...');
for xx=1:epoch
waitbar(xx/40, h, sprintf('已完成 0.03', i));
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
Diag3(1,:,xx)=mid(1:15);
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
    paraRQAPerio3(j,1,xx)=DET;
    paraRQAPerio3(j,2,xx)=ENTR;
    paraRQAPerio3(j,3,xx)=L;
    paraRQAPerio3(j,4,xx)=DETV;
    paraRQAPerio3(j,5,xx)=ENTRV;
end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j));
    paraCNPerio3(j,1,xx)=MATCHINGINDEX;
    paraCNPerio3(j,2,xx)=TRANSIVITY;
    paraCNPerio3(j,3,xx)=GLOBALEFFICIENCY;
    paraCNPerio3(j,4,xx)=MEANCLUSTER;
    paraCNPerio3(j,5,xx)=0;

end

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTPerio3(j,:,xx)=f;
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
Diag3(2,:,xx)=mid(1:15);
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
    paraRQAQuasi3(j,1,xx)=DET;
    paraRQAQuasi3(j,2,xx)=ENTR;
    paraRQAQuasi3(j,3,xx)=L;
    paraRQAQuasi3(j,4,xx)=DETV;
    paraRQAQuasi3(j,5,xx)=ENTRV;

end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j));
    paraCNQuasi3(j,1,xx)=MATCHINGINDEX;
    paraCNQuasi3(j,2,xx)=TRANSIVITY;
    paraCNQuasi3(j,3,xx)=GLOBALEFFICIENCY;
    paraCNQuasi3(j,4,xx)=MEANCLUSTER;
    paraCNQuasi3(j,5,xx)=0;
end

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTQuasi3(j,:,xx)=f;
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
Diag3(3,:,xx)=mid(1:15);
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
    paraRQAChao3(j,1,xx)=DET;
    paraRQAChao3(j,2,xx)=ENTR;
    paraRQAChao3(j,3,xx)=L;
    paraRQAChao3(j,4,xx)=DETV;
    paraRQAChao3(j,5,xx)=ENTRV;

end

for j=1:scale
    [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j).*Birecur(:,:,j)');
    paraCNChao3(j,1,xx)=MATCHINGINDEX;
    paraCNChao3(j,2,xx)=TRANSIVITY;
    paraCNChao3(j,3,xx)=GLOBALEFFICIENCY;
    paraCNChao3(j,4,xx)=MEANCLUSTER;
    paraCNChao3(j,5,xx)=0;

end

for j=1:scale
    [Dist] = RTA(Birecur(:,:,j));
    [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
    paraRTChao3(j,:,xx)=f;
end

xx

end
%% Random 
paraRQARandom=zeros(scale,5,epoch);
paraCNRandom=zeros(scale,5,epoch);
paraRTRandom=zeros(scale,2*Kn,epoch);
DiagRandom=zeros(1,40,epoch);

RR=0.05;
k=1;
h = waitbar(0, '开始处理...');
for k=1:epoch
    waitbar(i/40, h, sprintf('已完成0.05', i));
    y=randn(1,Size);
    distM=zeros(Size);
    for i=1:Size
        for j=1:Size
            distM(i,j)=abs(y(i)-y(j));
        end
    end
    [u,s,v]=svd(distM);
    mid=diag(s);
    DiagRandom(1,:,k)=mid(1:40);
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
        paraRQARandom(j,1,k)=DET;
        paraRQARandom(j,2,k)=ENTR;
        paraRQARandom(j,3,k)=L;
        paraRQARandom(j,4,k)=DETV;
        paraRQARandom(j,5,k)=ENTRV;

    end

    for j=1:scale
        [MATCHINGINDEX,MEANDEGREE,TRANSIVITY,GLOBALEFFICIENCY,MEANCLUSTER]= CN(Birecur(:,:,j).*Birecur(:,:,j)');
        paraCNRandom(j,1,k)=MATCHINGINDEX;
        paraCNRandom(j,2,k)=TRANSIVITY;
        paraCNRandom(j,3,k)=GLOBALEFFICIENCY;
        paraCNRandom(j,4,k)=MEANCLUSTER;
        paraCNRandom(j,5,k)=0;

    end

    for j=1:scale
        [Dist] = RTA(Birecur(:,:,j));
        [f,xi]=ksdensity(Dist,[0:0.5:Kn-0.5]);
        paraRTRandom(j,:,k)=f;
    end
end


