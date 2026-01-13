%% plot the S
figure(1)
subplot(2,3,1)
singularPerio=Diag(1,:,:);
meanSingularPerio=mean(singularPerio,3);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
fill([1:15,flip([1:15])],log10([meanSingularPerio+standSingularPerio,flip(meanSingularPerio-standSingularPerio)])...
    ,[0.5,0.5,0.5],'LineStyle','none','FaceAlpha',0.5)
hold on
plot([1:15],log10(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],log10(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('周期序列的奇异值分布图')

subplot(2,3,2)
singularPerio=Diag(2,:,:);
meanSingularPerio=mean(singularPerio,3);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
fill([1:15,flip([1:15])],log10([meanSingularPerio+standSingularPerio,flip(meanSingularPerio-standSingularPerio)])...
    ,[0.5,0.5,0.5],'LineStyle','none','FaceAlpha',0.5)
hold on
plot([1:15],log10(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],log10(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('准周期序列的奇异值分布图')

subplot(2,3,3)
singularPerio=Diag(3,:,:);
meanSingularPerio=mean(singularPerio,3);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
fill([1:15,flip([1:15])],log10([meanSingularPerio+standSingularPerio,flip(meanSingularPerio-standSingularPerio)])...
    ,[0.5,0.5,0.5],'LineStyle','none','FaceAlpha',0.5)
hold on
plot([1:15],log10(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],log10(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('混沌序列的奇异值分布图')

subplot(2,3,4)
singularPerio=Diag(1,:,:);
meanSingularPerio=mean(singularPerio,3);
meanSingularPerio=cumsum(meanSingularPerio)/sum(meanSingularPerio);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
periocumsum=meanSingularPerio;
plot([1:15],(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('周期序列的奇异值分布图')

subplot(2,3,5)
singularPerio=Diag(2,:,:);
meanSingularPerio=mean(singularPerio,3);
meanSingularPerio=cumsum(meanSingularPerio)/sum(meanSingularPerio);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
quasiperiocumsum=meanSingularPerio;
plot([1:15],(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('准周期序列的奇异值分布图')

subplot(2,3,6)
singularPerio=Diag(3,:,:);
meanSingularPerio=mean(singularPerio,3);
meanSingularPerio=cumsum(meanSingularPerio)/sum(meanSingularPerio);
varSingularPerio=var(singularPerio,1,3);
standSingularPerio=sqrt(varSingularPerio);
chaocumsum=meanSingularPerio;
plot([1:15],(meanSingularPerio),'Color',[0,0,0],'LineStyle','--','LineWidth',1)
hold on
scatter([1:15],(meanSingularPerio),'black','o')
grid on
xlabel('回复阶数')
ylabel('log 10 奇异值')
title('混沌序列的奇异值分布图')
%% analyse the recover under R=0.05
load('E:\研究生\repo\Rrecuurence-matrix-eigenface\data\approach1\R0.05.mat')
figure(2)
subplot(3,2,1)
para=mean(paraRQAPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','DETV','ENTRV')
title('RQA PErio')
subplot(3,2,2)
para=mean(paraCNPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN PErio')
subplot(3,2,3)
para=mean(paraRQAQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','DETV','ENTRV')
title('RQA quasiPErio')
subplot(3,2,4)
para=mean(paraCNQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN quasiperio')
subplot(3,2,5)
para=mean(paraRQAChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','DETV','ENTRV')
title('RQA chao')
subplot(3,2,6)
para=mean(paraCNChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN chao')
%% analyse the recover under R=0.1
load('E:\研究生\repo\Rrecuurence-matrix-eigenface\data\approach1\R0.1.mat')
figure(3)
subplot(3,2,1)
para=mean(paraRQAPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
% plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
% hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','ENTRV')
title('RQA PErio')
subplot(3,2,2)
para=mean(paraCNPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN PErio')
subplot(3,2,3)
para=mean(paraRQAQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
% plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
% hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','ENTRV')
title('RQA quasiPErio')
subplot(3,2,4)
para=mean(paraCNQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN quasiperio')
subplot(3,2,5)
para=mean(paraRQAChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','DETV','ENTRV')
title('RQA chao')
subplot(3,2,6)
para=mean(paraCNChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot([1:12],relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot([1:12],relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot([1:12],relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN chao')
%%  lineral regression under RR=0.1
res=1-[chaocumsum(1:12);chaocumsum(1:12);quasiperiocumsum(1:12)...
    ;quasiperiocumsum(1:12);periocumsum(1:12);periocumsum(1:12)];
figure(4)
subplot(3,2,1)
para=mean(paraRQAPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(6,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
% plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
% hold on
plot(res(6,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','ENTRV')
title('RQA PErio')
subplot(3,2,2)
para=mean(paraCNPerio,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(6,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot(res(6,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN PErio')
subplot(3,2,3)
para=mean(paraRQAQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(4,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
% plot([1:12],relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
% hold on
plot(res(4,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','ENTRV')
title('RQA quasiPErio')
subplot(3,2,4)
para=mean(paraCNQuasi,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(4,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot(res(4,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN quasiperio')
subplot(3,2,5)
para=mean(paraRQAChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(2,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('DET','ENTR','L','DETV','ENTRV')
title('RQA chao')
subplot(3,2,6)
para=mean(paraCNChao,3);
relaErr=zeros(12,5);
for i=1:12
    relaErr(i,:)=abs(para(i,:)-para(12,:))./para(12,:);
end
ref=relaErr(1,:);
for i=1:12
    relaErr(i,:)=relaErr(i,:)./ref;
end
plot(res(2,:),relaErr(:,1),'Color',[1,0,0],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,2),'Color',[0,1,0],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,3),'Color',[0,0,1],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,4),'Color',[0.5,0.5,0.5],'LineStyle','--')
hold on
plot(res(2,:),relaErr(:,5),'Color',[0.5,0.25,0.75],'LineStyle','--')
legend('MATCHINGINDEX','TRANSIVITY','GLOBALEFFICIENCY','MEANCLUSTER')
title('CN chao')


