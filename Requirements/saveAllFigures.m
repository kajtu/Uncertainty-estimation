function saveAllFigures(plotFolderName,saveFigs,savePDF)

if nargin < 2 || (nargin >1 && saveFigs)
    disp('Saving figures...')    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for fig = 1:length(FigList)
        currFig = FigList(fig);
        figName = currFig.Name;
        print(currFig,fullfile(plotFolderName,figName),'-dpng','-r600')
        if nargin >2 && savePDF
            exportgraphics(currFig,[fullfile(plotFolderName,figName) '.pdf'],'ContentType','vector')
        end
        close(currFig)
        fprintf('%0.1f \n',100*fig/length(FigList))
    end
else
    disp('Not saving any figures.')
end

end