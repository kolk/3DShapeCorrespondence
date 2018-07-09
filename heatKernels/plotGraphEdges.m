function plotGraphEdges(vertCoords, adjMatrix)
   n = size(adjMatrix, 2);
   figure;
   hold on;
   scatter(vertCoords(:,2),vertCoords(:,1), 'k*');
   for i = 1: n
       idx = find(adjMatrix(i,:)>0);
       for j = 1:length(idx)
           plot([vertCoords(idx(j),2), vertCoords(i,2)], [vertCoords(idx(j),1), vertCoords(i,1)], 'b-');
       end
   end
   hold off;
end