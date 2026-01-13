function Dist = RTA(RP)
%UNTITLED 此处显示有关此函数的摘要
%   recurrence time analyses represent the time of recurrence
%   analyse the mean RT under different lines
Size=size(RN,1);
Dist=zeros(Size,1);
for i=1:Size
    mid=RN(:,i);
    mid=[0;mid;0];
    Diff=diff(mid);
    place=find(Diff==1);
    place1=[0;place];
    place2=[place;0];
    distPlace=place2-place1;
    distPlace=distPlace(2:size(place1)-1);
    Dist(i)=mean(distPlace);
end
end