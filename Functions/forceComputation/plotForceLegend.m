function plotForceLegend
figure; 
set(gcf, 'Units', 'centimeter', 'OuterPosition', [1 1 10 7]);
hold on
quiverLegend = [4 3 -2 0; 5 3 2 0; 4 2 -2 0; 5 2 2 0]; % create the vectors
QuiverLegendColor = [0.5 0 0.1; 0.8 0.5 0.5; 0 0.2 0.5;0.45 0.68 0.8]; % load the corresponding colors

% Plot the colored vectors
for i = 1:size(quiverLegend,1)
    quiver(quiverLegend(i,1),quiverLegend(i,2),quiverLegend(i,3),quiverLegend(i,4),'MaxHeadSize',0.7,'LineWidth', 2,'Color', QuiverLegendColor(i,:)) % reference vector scale
end

% Add thrust and drag text
text(2.75,1.5,'Drag','Color',[0 0 0],'FontSize',18) % change the text and its color as needed
text(5.25,1.5,'Thrust','Color',[0 0 0],'FontSize',18) % change the text and its color as needed

% Add pressure contribution text
text(1.5,3,'(+)','Color',[0 0 0],'FontSize',18) % change the text and its color as needed
text(1.5,2,'(-)','Color',[0 0 0],'FontSize',18) % change the text and its color as needed

% Add pressure label
text(0.8,3.5,'Pressure','Color',[0 0 0],'FontSize',18) % change the text and its color as needed

axis equal
axis off
axis([0.5 7 1 4])
formatFigure
end