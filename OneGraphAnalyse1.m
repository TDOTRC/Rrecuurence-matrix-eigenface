%% cauculate the PCA based on one matrix
traitVec=zeros(Size,Vol*9);
%% lorentz stability under approach 2
load("E:\研究生\论文1\论文1\数据\lor距离\data4.mat")
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(data4(3*j,3*(i-1)+1)-data4(3*k,3*(i-1)+1));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(log(abs(norm(:,:,i))))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+0)=eVec(:,1,i);
end
%% sin random
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
t=[0.01:0.01:20];
z=zeros(Vol,Size);
for i=1:Vol
 f=rand(1);
 f=f+0.5;
z(i,:)=sin(2*pi*f*t);
end
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(z(i,k)-z(i,j));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(abs(norm(:,:,i)))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+Vol)=eVec(:,1,i);
end
%% rossler periodic
load("E:\研究生\论文1\论文1\数据\lor距离\data3.mat")
data3=x;
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
t=[0.01:0.01:20];
z=zeros(Vol,Size);
for i=1:Vol
 f=rand(1);
 f=f+0.5;
z(i,:)=sin(6*f*t);
end
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(z(i,k)-z(i,j));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(abs(norm(:,:,i)))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+2*Vol)=eVec(:,1,i);
end
%%  gaussed lorentz stability under approach 2
load("E:\研究生\论文1\论文1\数据\lor距离\data4.mat")
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
Rand=randn(Size,1);
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(data4(3*j,3*(i-1)+1)-data4(3*k,3*(i-1)+1)+Rand(j)-Rand(k));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(log(abs(norm(:,:,i))))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+3*Vol)=eVec(:,1,i);
end
%% gaussed sin
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
t=[0.01:0.01:20];
z=zeros(Vol,Size);
Rand=randn(Size,1);
for i=1:Vol
 f=rand(1);
 f=f+0.5;
z(i,:)=sin(2*pi*f*t);
end
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(z(i,k)+Rand(j)-z(i,j)-Rand(k));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(abs(norm(:,:,i)))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+4*Vol)=eVec(:,1,i);
end
%% gaussed rossler
load("E:\研究生\论文1\论文1\数据\lor距离\data3.mat")
data3=x;
Size=2000;
Col=15;
Vol=45;
distance=zeros(Size,Size,Vol);
pca=zeros(Size,Col,Vol);
t=[0.01:0.01:20];
z=zeros(Vol,Size);
Rand=randn(Size,1);
for i=1:Vol
 f=rand(1);
 f=f+0.5;
z(i,:)=sin(6*f*t);
end
for i=1:Vol
    for j=1:Size
        for k=1:Size
            distance(j,k,i)=abs(z(i,k)+Rand(j)-Rand(k)-z(i,j));
        end
    end
    i
end
Cov=zeros(Size,Size);
eVec=zeros(Size,Col,Vol);
eVal=zeros(Col,Vol);
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
    mid=flip(diag(eValue));
    eVal(:,i)=mid(1:15)';
    i
end
%verify thewhether it is normal
norm=zeros(Col,Col,Vol);
for i=1:Vol
        norm(:,:,i)=eVec(:,:,i)'*eVec(:,:,5);
end
figure(1)
for i=1:10
    subplot(2,5,i)
    imagesc(abs(norm(:,:,i)))
end
figure(2)
for j=1:10
    subplot(2,5,j)
    scatter([1:15],log(eVal(:,j)),4,"red","filled")
end
for i=1:Vol
    traitVec(:,i+5*Vol)=eVec(:,1,i);
end