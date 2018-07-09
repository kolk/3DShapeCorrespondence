function drawCorrespondences(data, label, labels, correspondences, k, numShapes, numNodes, cliqueWeight)
    map = [1, 1, 0
    0, 1, 1
    1, 0, 1
    1, 0, 0
    0, 1, 0
    0, 0, 1
    0.5, 0.5, 0.5
    0.2, 0.2, 0.2
    0.1, 0.8, 0.4
    1, 1, 1
    0.8 0.5 0.5
    0.3 0.7 0.9
    0.6 0.6 1
    0.3 0.2 0.7
    0.7 0.1 0.1
    ];
    newLabel = zeros(numNodes(1), numShapes);
    for gauss = 1:k
        shape = 1;
        rw = find(correspondences(:,1) == gauss);
        tLabel = find(label(:,shape) == gauss);
        lbl = gauss;%correspondences(gauss,1);
        newLabel(tLabel,shape) = lbl;


        shape = 2;
        tLabel = find(label(:,shape) == correspondences(rw,2));
        newLabel(tLabel, shape) = lbl;

        if numShapes > 2
            shape = 3;
            %find shape 2 label in 1st column
            correspondences(rw,2)+(shape-2)*k
            corrIndx = find(correspondences(:,1) == correspondences(rw,2)+(shape-2)*k)
            gaussNum = correspondences(corrIndx,2)-(shape-2)*k
            tLabel = find(label(:, shape) == gaussNum);
            newLabel(tLabel, shape) = lbl;

            if numShapes > 3
                shape = 4;
                corrIndx = find(correspondences(:,1) == gaussNum+(shape-2)*k);
                gaussNum = correspondences(corrIndx,2)-(shape-2)*k;
                tLabel = find(label(:, shape) == gaussNum);
                newLabel(tLabel, shape) = lbl;
            end
        end
    end
    %end


    figure;
    colormap(map(1:k,:));
    trisurf(data(1).surface.TRIV, data(1).surface.X, data(1).surface.Y, data(1).surface.Z, newLabel(:,1));
    title(['Weight of Clique ' num2str(cliqueWeight)]);
    %title(['cost of clique ' num2str(affinityofClique)]);
    %lcolorbar(labels(1:k));
    translate = 150;
    for i = 2:numShapes
        hold on;
        disp(['displaying surface ' num2str(i)]);
        trisurf(data(i).surface.TRIV, data(i).surface.X + repmat([translate], size(data(i).surface.X, 1), 1), data(i).surface.Y, data(i).surface.Z, newLabel(:,i));
        translate = translate + 100;
    end
end