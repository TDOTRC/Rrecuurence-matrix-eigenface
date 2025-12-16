%% use kmeans approach to analyse the classification
%% pre cook the data
ortho=cell2mat(struct2cell(load("matlab.mat")));
parad=cell2mat(struct2cell(load('matlabTiled.mat')));
ortho=ortho';
parad=parad';
%% orthodox condition
% set the super para
k=9;
epoch=100;
repli=5;
[idx, C, sumd] = kmeans(ortho, k,'MaxIter', epoch,'Replicates', repli);
plot(C(9,:))
%% paradox condition
k=9;
epoch=100;
repli=5;
[idx, C, sumd] = kmeans(parad, k,'MaxIter', epoch,'Replicates', repli);
plot(C(9,:))