K = 2;
pic=im2double(imread('ponyo.jpeg'));
Reshaped = reshape(pic,size(pic,2)*size(pic,1),3);
D = zeros(size(Reshaped,1),K+2);
C = Reshaped(ceil(rand(K,1)*size(Reshaped,1)),:);
I = zeros(size(Reshaped));

for i = 1:15
    for j = 1:size(Reshaped,1)
        for n = 1:K
            D(j,n) = norm(Reshaped(j,:) - C(n,:));
        end
        [D(j,K+2),D(j,K+1)] = min(D(j,1:K));
    end
    
    for j = 1:K
        C(j,:) = mean(Reshaped(D(:,K+1) == j,:));
        if sum(isnan(C(:))) ~=0
            for l = 1:size(find(isnan(C(:,1)) == true),1)
                C(N(l),:) = Reshaped(randi(size(Reshaped,1)),:);
            end
        end
    end
end

for i = 1:K
    index = find(D(:,K+1)==i);
    I(index,:) = repmat(C(i,:),size(index,1),1);
end

result = reshape(I,size(pic,1),size(pic,2),3);
figure
imshow(result)
title(['K=',num2str(K)]);