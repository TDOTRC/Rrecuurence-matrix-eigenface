%% cauculate the PCA based on one matrix
%% sin
t=[0.01:0.01:20];
y=0.1*sin(2*pi*t);
z=y;
distance=zeros(2000,2000);
for i=1:2000
    for j=1:2000
        distance(i,j)=abs(z(i)-z(j));
    end
end
figure(2)
deposit=zeros(20,2000);
for k=1:20
    deposit(k,:)=distance(5*k,:);
end
clear k
for k=1:20
    plot(t,deposit(k,:))
    hold on
end
ave=zeros(2000,1);
ave=sum(distance,2)/2000;
for i=1:2000
    distance(:,i)=distance(:,i)-ave;
end
average=mean(distance,2);
cov=(distance)*(distance)';
[u,s,v]=svd(distance,'econ');
[eVec,eValue]=eig(cov);
mid=eVec(:,1996);
plot(mid);
ans=fft(mid);
plot([1:2000],ans)



