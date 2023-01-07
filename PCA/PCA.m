clc;
load('1d_pcadata.mat');
Mean = mean(X,1);
Y = bsxfun(@minus,X,Mean);
[coeff,score,latent] = pca(Y);
[V,D] = eig(cov(Y));
%
figure
scatter(X(:,1), X(:,2), 25,'LineWidth',1.5)
xlim([0 7]); ylim([2 8])
hold on;
quiver(Mean(1),Mean(2),V(1,2),V(2,1),latent(1),'k','LineWidth',1.2);
quiver(Mean(1),Mean(2),V(2,2),V(1,1),latent(2),'k','LineWidth',1.2);
hold off;
title('PCA')
%}

load('faces.mat')
MeanFace = mean(X,1);
Y = bsxfun(@minus,X,MeanFace);
[coeff,score,latent] = pca(Y);
[V,D] = eig(cov(Y),'vector');
[D, x] = sort(D, 'descend');
D(1:5)
V = V(:, x);
Eigen = V(:, 1:36);
Eigen = reshape(Eigen,[32,32,36]);
figure
montage(Eigen, 'DisplayRange', [])
title('Eigenfaces')

D = V';
AllFace = reshape(D' * D * Y' + repmat(MeanFace(:), 1, 5000), [32, 32, 5000]);
D = V(:,1:20)';
TwentyFace = reshape(D' * D * Y' + repmat(MeanFace(:), 1, 5000), [32, 32, 5000]);

figure; montage(AllFace(:,:,1:100), 'DisplayRange', [])
title('All Eigenvectors')
figure; montage(TwentyFace(:,:,1:100), 'DisplayRange', [])
title('20 Eigenvectors')
%}