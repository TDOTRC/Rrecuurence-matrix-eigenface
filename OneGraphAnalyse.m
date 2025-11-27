%% cauculate the PCA based on one matrix
%% lorentz stability under approach 2
Size=2000;
Col=15;
Vol=20;
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
for i=1:20
    Cov=distance(:,:,i)*distance(:,:,i)';
    [eVector,eValue]=eig(Cov);
    for j=1:Col
        eVec(:,j,i)=eVector(:,Size-(j-1));
    end
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