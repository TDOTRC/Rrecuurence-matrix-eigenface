%These codes are for deciding which kinds of"faces" shall be 
% included in to the data set.
%% Parameters predecision
SizeDM=2000;
t_span=200;
dataSet=zeros(SizeDM^2,840);
%% cauculate the lorentz based graph 1-70
D=zeros(SizeDM,SizeDM);
sigma=10;
beta=2.66;
for i=1:10
    rho=28+i;
    for j=1:7
    initial_point=rand(1,3);
    [t, trajectory] = lorenz_trajectory(sigma,rho,beta,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1);
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+0)=D(:);
    end
    i
end
%% cauculate the rossler based graph 71-140
D=zeros(SizeDM,SizeDM);
b=0.2;
c=5.7;
for i=1:10
    a=0.2+0.01*i;
    for j=1:7
    initial_point=rand(1,3);
    [t, trajectory] = rossler_trajectory(a,b,c,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1);
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+1*70)=D(:);
    end
    i
end
%% cauculate the sin based graph 141-210
D=zeros(SizeDM,SizeDM);
for i=1:10
    A=0.1*i;
    for j=1:7
        f=0.1*j;
        timeSeries=A*sin(f*[0.1:0.1:0.1*SizeDM]);
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+2*70)=D(:);
    end
    i

end
%% cauculate the sined lorentz based graph 211-280
D=zeros(SizeDM,SizeDM);
sigma=10;
beta=2.66;
rho=28;
for i=1:10
    sined=5*sin(2*pi*0.001*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    initial_point=rand(1,3);
    [t, trajectory] = lorenz_trajectory(sigma,rho,beta,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1)+sined';
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+3*70)=D(:);
    end
    i
end
%% cauculate the sined rossler based graph 281-350
D=zeros(SizeDM,SizeDM);
a=0.2;
b=0.2;
c=5.7;
for i=1:10
    sined=5*sin(2*pi*0.001*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    initial_point=rand(1,3);
    [t, trajectory] = rossler_trajectory(a,b,c,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1)+sined';
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+4*70)=D(:);
    end
    i
end
%% cauculate the sined sin graph 351-420
D=zeros(SizeDM,SizeDM);
for i=1:10
    sined=5*sin(2*pi*0.01*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    f=0.1*j;
    timeSeries=15*sin([0.1:0.1:0.1*SizeDM])+sined;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+5*70)=D(:);
    end
    i
end
%% cauculate the drifted lorentz graph 421-490
D=zeros(SizeDM,SizeDM);
sigma=10;
beta=2.66;
for i=1:10
    rho=28+i;
    for j=1:7
    drift=0.01*j*[1:1:2000]';
    initial_point=rand(1,3);
    [t, trajectory] = lorenz_trajectory(sigma,rho,beta,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1);
    timeSeries=timeSeries+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+6*70)=D(:);
    end
    i
end
%% cauculate the drifted rossler graph 491-560
D=zeros(SizeDM,SizeDM);
a=0.2;
b=0.2;
c=5.7;
for i=1:10
    for j=1:7
    initial_point=rand(1,3);
    drift=0.005*j*[1:1:2000]';
    [t, trajectory] = rossler_trajectory(a,b,c,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1)+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+7*70)=D(:);
    end
    i
end
%% cauculate the drifted sin graph 561-630
D=zeros(SizeDM,SizeDM);
for i=1:10
   A=4;
    for j=1:7
        f=0.1*j;
        drift=0.01*j*[1:1:2000];
        timeSeries=A*sin(f*[0.1:0.1:0.1*SizeDM])+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+8*70)=D(:);
    end
    i

end
%% cauculate the drifted sined lorentz graph 631-700
D=zeros(SizeDM,SizeDM);
sigma=10;
beta=2.66;
rho=28;
for i=1:10
    sined=5*sin(2*pi*0.001*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    drift=0.01*j*[1:1:2000]';
    initial_point=rand(1,3);
    [t, trajectory] = lorenz_trajectory(sigma,rho,beta,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1)+sined'+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+9*70)=D(:);
    end
    i
end
%% cauculate the drifted sined rossler graph 701-770
D=zeros(SizeDM,SizeDM);
a=0.2;
b=0.2;
c=5.7;
for i=1:10
    sined=5*sin(2*pi*0.001*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    drift=0.005*j*[1:1:2000]';
    initial_point=rand(1,3);
    [t, trajectory] = rossler_trajectory(a,b,c,initial_point, t_span);
    timeSeries=trajectory([1:8:16000],1)+sined'+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+10*70)=D(:);
    end
    i
end
%% cauculate the drifted sined sin graph 771-840
D=zeros(SizeDM,SizeDM);
for i=1:10
    sined=5*sin(2*pi*0.001*(1+0.1*i)*[1:SizeDM]);
    for j=1:7
    f=0.1*j;
    drift=0.01*j*[1:1:2000];
    timeSeries=2.5*sin([0.1:0.1:0.1*SizeDM])+sined+drift;
    for k=1:SizeDM
        for m=1:SizeDM
            D(k,m)=abs(timeSeries(k)-timeSeries(m));
        end
    end
    dataSet(:,i*j+11*70)=D(:);
    end
    i
end
%% check whether dataset is real

%% pre cook the dataset(recuce the average)
dataMean=zeros(1,SizeDM^2);
for m=1:SizeDM^2
    dataMean(m)=mean(dataSet(m,:));
    dataSet(m,:)=dataSet(m,:)-mean(dataSet(m,:))*ones(1,840);
   
end
dataMean=dataMean';
dataAve=reshape(dataMean,SizeDM,SizeDM);
dataCooked=log(dataAve+diag(ones(SizeDM)));
dataCooked=dataCooked/max(dataCooked,[],'all');
figure(1)
imshow(dataCooked)
%% determine the cov matrix and diconpose it
V=dataSet'*dataSet;
[enginVec,enginValue]=eig(V);
verify=enginVec*enginVec';
mid=diag(enginValue);
tolerant=1*10^(-5);
for i=1:840
    m=mid(i);
    if m<tolerant
        mid(i)=0;
    else
        mid(i)=mid(i);
    end
end
enginValue=flip(mid);
engin=log10(enginValue+eps*ones(1,840));
middle=cumsum(enginValue)/sum(enginValue);
figure(2)
subplot(1,2,1)
plot(engin)
xlabel('特征值序号')
ylabel('log特征值')
title('特征值分布图')
subplot(1,2,2)
scatter([1:100],middle(1:100),'black')
xlabel('特征值序号')
ylabel('相对恢复占比')
hold on
plot([1:100],0.99,'r')
title('特征值恢复状况图')
%% cauculate the eigen face and relavent varience
Consider=90;
relaventEvalue=zeros(1,Consider);
eigenFace=zeros(SizeDM^2,Consider);
for i=0:Consider-1
    relaventEvalue(i+1)=mid(840-i);
    eigenFace(:,i+1)=dataSet*enginVec(:,840-i);
end
figure(3)
subplot(3,2,1)
imagesc(reshape(eigenFace(:,1),SizeDM,SizeDM))
title('1阶engin RP')
subplot(3,2,2)
imagesc(reshape(eigenFace(:,2),SizeDM,SizeDM))
title('2阶engin RP')
subplot(3,2,3)
imagesc(reshape(eigenFace(:,3),SizeDM,SizeDM))
title('3阶engin RP')
subplot(3,2,4)
imagesc(reshape(eigenFace(:,4),SizeDM,SizeDM))
title('4阶engin RP')
subplot(3,2,5)
imagesc(reshape(eigenFace(:,5),SizeDM,SizeDM))
title('5阶engin RP')
subplot(3,2,6)
imagesc(reshape(eigenFace(:,6),SizeDM,SizeDM))
title('6阶engin RP')

figure(8)
subplot(3,2,1)
imagesc(reshape(eigenFace(:,7),SizeDM,SizeDM))
title('7阶engin RP')
subplot(3,2,2)
imagesc(reshape(eigenFace(:,8),SizeDM,SizeDM))
title('8阶engin RP')
subplot(3,2,3)
imagesc(reshape(eigenFace(:,9),SizeDM,SizeDM))
title('9阶engin RP')
subplot(3,2,4)
imagesc(reshape(eigenFace(:,10),SizeDM,SizeDM))
title('10阶engin RP')
subplot(3,2,5)
imagesc(reshape(eigenFace(:,11),SizeDM,SizeDM))
title('11阶engin RP')
subplot(3,2,6)
imagesc(reshape(eigenFace(:,12),SizeDM,SizeDM))
title('12阶engin RP')
%cauculate the frequent domain
figure(4)
subplot(2,3,1)
a1=fft(eigenFace(1:SizeDM,1));
plot(log10(abs(a1(1:0.25*SizeDM))))
title('1阶engin RP fft')
grid on
subplot(2,3,2)
title('2阶engin RP fft')
a1=fft(eigenFace(1:SizeDM,2));
plot(log10(abs(a1(1:0.25*SizeDM))))
grid on
subplot(2,3,3)
a1=fft(eigenFace(1:SizeDM,3));
plot(log10(abs(a1(1:0.25*SizeDM))))
title('3阶engin RP fft')
grid on
subplot(2,3,4)
a1=fft(eigenFace(1:SizeDM,4));
plot(log10(abs(a1(1:0.25*SizeDM))))
title('4阶engin RP fft')
grid on
subplot(2,3,5)
a1=fft(eigenFace(1:SizeDM,5));
plot(log10(abs(a1(1:0.25*SizeDM))))
title('5阶engin RP fft')
grid on
subplot(2,3,6)
a1=fft(eigenFace(1:SizeDM,6));
plot(log10(abs(a1(1:0.25*SizeDM))))
title('6阶engin RP fft')
grid on

figure(5)
subplot(2,3,1)
plot(eigenFace(2*SizeDM+1:3*SizeDM,1))
title('1阶engin RP')
subplot(2,3,2)
plot(eigenFace(2*SizeDM+1:3*SizeDM,2))
title('2阶engin RP')
subplot(2,3,3)
plot(eigenFace(2*SizeDM+1:3*SizeDM,3))
title('3阶engin RP')
subplot(2,3,4)
plot(eigenFace(2*SizeDM+1:3*SizeDM,4))
title('4阶engin RP')
subplot(2,3,5)
plot(eigenFace(2*SizeDM+1:3*SizeDM,5))
title('5阶engin RP')
subplot(2,3,6)
plot(eigenFace(2*SizeDM+1:3*SizeDM,6))
title('6阶engin RP')
% a little try
finaldata=dataSet(:,113);
position=zeros(1,Consider);
for i=1:Consider
    position(i)=eigenFace(:,i)'*finaldata;
end
composeData=zeros(SizeDM^2,1);
for j=1:Consider
composeData=composeData+position(j)*eigenFace(:,j);
end
 figure(6)
 subplot(1,2,1)
 ans1=reshape(finaldata,SizeDM,SizeDM);
 imagesc(ans1)
 subplot(1,2,2)
 ans2=(reshape(composeData,SizeDM,SizeDM))+eps*diag(SizeDM);
 ansFinal=ans2*10^-11;
 imagesc(ansFinal)
%normal
 d=zeros(Consider,Consider);
 for i=1:Consider
     for j=1:Consider
         d(i,j)=eigenFace(:,i)'*eigenFace(:,j);
     end
 end
 figure(7)
 imagesc(log(abs(d)))
%% verify the correct of distance matrix
