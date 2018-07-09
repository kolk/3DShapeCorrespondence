function [M, l, Dist] = associationGraphSequence3(mu, sigma, numSequence, k, covthresh)
        % k = number of features, d = dimensionality of data
        % model1.mu - k x d matrix; 
        % model1.Sigma - d x d x k covariance matrix
        % covthresh - threshold for eliminating the assignments of
        % assocation graph

         % size of all possible assignments
        n = k^numSequence;

        % list of candidate assignments
        l = zeros(n,numSequence);
        t = 1;
        
        covDiff = zeros(n,1);
        for i = 1:k
            for j = 1:k
                % check if variance difference threshold is met
                cov1 = sigma(:,:,i,1);
                cov2 = sigma(:,:,j,2);
                diff = (cov1 - cov2).^2;
                val1 = sum(diff(:));
                covDiff(t) = val1;
                l(t,1) = i;
                l(t,2) = j;
                t = t + 1;
            end
        end

        covDiff = (covDiff - min(covDiff))/(max(covDiff) - min(covDiff));
        [~, ind] = sort(covDiff);
        %disp(sprintf('%d %d %f\n',l(ind,1),l(ind,2),covDiff(ind)));
      
        %covthresh = mean(covDiff); 
        indx = find(covDiff <= covthresh);
        %indx = 1:int8(0.60*length(l));
        l = l(indx, :);
       
        
        n = length(l);
        %distance matrix
        Dist = zeros(n,n);
        covDist = zeros(n,n);
        size(Dist);
        
        % affinity matrix
        M = zeros(n, n);

        
        for a = 1:n
            for b = 1:a
                %d = -1;
                % measure of perservation of relative pairwise geometry
                if a ~= b
                    
                    dist = zeros(numSequence,1);
                    % find distance between 2 gaussian centers in shape i
                    for i = 1:numSequence
                        p = l(a,i);
                        q = l(b,i);
                        % mu(p,:,i) - mean of gaussian p of shape i
                        dist(i) = sqrt(sum((mu(p,:,i) - mu(q,:,i)).^2));
                    end
                       
                  
                    d = 0;
                   for i = 1:numSequence -1
                       d = d + abs(dist(i+1) - dist(i));
                   end
                   
                     %find forb norm of covariances
                   %for i = 1:numSequence-1
                        p = l(a,:);
                        q = l(b,:);
                        covDiff1 = (sigma(:,:,p(1),1) - sigma(:,:,p(2),2)).^2;
                        covDiff2 = (sigma(:,:,q(1),1) - sigma(:,:,q(2),2)).^2;
                        covDist1 = sqrt(sum(covDiff1(:)));
                        covDist2 = sqrt(sum(covDiff2(:)));
                   %end
                   
                   covDist(a,b) = covDist1 + covDist2;
                   covDist(b,a) = covDist1 + covDist2;
                   Dist(a,b) =  d; %+ 1*covDist1 + 1*covDist2;
                   Dist(b,a) = d; %+ 1*covDist1 + 1*covDist2;                
                else
                    %{
                    p = l(a,:);
                    covDiff1 = (sigma(:,:,p(1),1) - sigma(:,:,p(2),2)).^2;
                    covDist1 = sqrt(sum(covDiff1(:)));
                    covDist(a,b) = covDist1;
                    covDist(b,a) = covDist1;
                    Dist(a,b) = 0; %+ 1*covDist1 + 1*covDist2;
                    Dist(b,a) = 0; %+ 1*covDist1 + 1*covDist2; 
                    %}
                end

            end
        end
        
        %{
        for i = 1:size(Dist, 2)
           Dist(:,i) = (Dist(:,i)-min(Dist(:,i))) / (max(Dist(:,i)) - min(Dist(:,i)));
        end
        for i = 1:size(covDist, 2)
           covDist(:,i) = (covDist(:,i)-min(covDist(:,i))) / (max(covDist(:,i)) - min(covDist(:,i)));
        end
        %}
        
        %Dist = (Dist - min(Dist(:)))/(max(Dist(:)) - min(Dist(:)));
        %covDist = (covDist - min(covDist(:)))/(max(covDist(:)) - min(covDist(:)));
        
        M = zeros(n,n);
        alpha = 0.5;
        for i = 1:n
            for j = 1:n
                if Dist(i,j) > 0
                    M(i,j) = alpha*exp(-Dist(i,j)/(mean(mean(Dist))));
                end
                if covDist(i,j) > 0
                    M(i,j) = M(i,j) + (1-alpha)*exp(-covDist(i,j)/(mean(mean(covDist))));
                end                                
            end
        end
        
        %Dist =  0.5 * Dist + 0.5 * covDist;
        
        %{
        for p = 1: size(Dist,1)
            Dist(p,:) = Dist(p,:) ./ sum(Dist(p,:));
        end
        %}
        %Dist = (Dist - min(Dist(:)))/(max(Dist(:)) - min(Dist(:)));
%         M = zeros(n,n);
%         sortedDist = sort(Dist);
%         sortedDist(1:10)
%         for i = 1:n
%             for j = 1:n
%                 if i ~= j
%                     if(Dist(i,j) <= 0.01)
%                         M(i,j) = 1 / 0.01;
%                     
%                     else
%                         M(i,j) = 1/Dist(i,j);
%                     end
%                 end
%              end
%         end
% 
       % normalize the matrix
       
       %{
       for i = 1:size(M, 2)
           M(i,:) = %(M(i,:)-min(M(i,:))) / (max(M(i,:)) - min(M(i,:)));
       end
         %}
       
%      figure;
%      imagesc(Dist);
%      title('Dist');
       
end
