function [pval,bestcost,fullParamVector] = extractPLparamfromfolder(folderPath,numParams,paramNamesALL)
% For a PL folder, find the parameter set with lowest cost and extract
% which parameter value the fixed parameter had, which the best cost was,
% and the full parameter vector including the fixed parameter.

        if ~isempty(dir([folderPath,'/opt-*']))
            [param,bestcost] = findBestParams(folderPath,0,numParams-1,0);
            if ~isnan(param)
                %add the parameter that was fixed
                folderName = split(folderPath,'/');
                folderName = folderName{end};
                folderName = split(folderName,'\');
                folderName = folderName{end};
                pnamefull = split(folderName,'_');
                %parameter name
                pname = strcat(pnamefull(1:end-1),'_');
                pname = strcat(pname{:});
                pname = pname(1:end-1);
                pind = find(strcmp(paramNamesALL,pname));
                %parameter value
                pval = pnamefull{end};
                pval = str2double(pval);
                fullParamVector = zeros(1,length(param)+1);
                fullParamVector(pind) = pval;
                fullParamVector(1:pind-1) = param(1:pind-1);
                fullParamVector(pind+1:end) = param(pind:end);
                if numParams == length(fullParamVector)+4 %add corr parameters
                    fullParamVector = [fullParamVector 40 40 40 40];
                end
            else
                pval=NaN;bestcost=NaN;fullParamVector=nan(1,numParams);
            end
        else
            pval=NaN;bestcost=NaN;fullParamVector=nan(1,numParams);
        end

end