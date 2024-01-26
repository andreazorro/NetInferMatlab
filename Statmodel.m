function [sNet_st] = Statmodel(gene_names, regulators, expressiondata)

    distribution = {'Beta','Binomial','BirnbaumSaunders','Burr','Exponential',...
        'ExtremeValue','Gamma','GeneralizedExtremeValue','GeneralizedPareto',...
        'HalfNormal','InverseGaussian','Kernel','Logistic','Loglogistic',...
        'Lognormal','Nakagami','NegativeBinomial','Normal','Poisson','Rayleigh',...
        'Rician','Stable','tLocationScale','Weibull'};
    
    tfs = find(ismember(gene_names,regulators));
    
    tfexpression = expressiondata(:,tfs);
    
    ngenes = size(expressiondata,2);
    ntf = size(tfs,2);
    ndist = size(distribution,2);
    
    expressiondata2 = distributed(expressiondata);
    
    spmd
        
        localdata = getLocalPart(expressiondata2);

        pv = zeros(ntf,size(localdata,2));
        
        for j = 1:size(localdata,2)
    
            Y = localdata(:,j);
            meanY = mean(Y);
            maxY = max(Y);
            minY = min(Y);
            stY = (Y - meanY)/(maxY - minY);
            stYm = stY+abs(min(stY))+1;
    
            h = zeros(ndist,1);
            p = zeros(ndist,1);

            for k = 1:ndist
                warning ('off','all')
                try
                    pd = fitdist(stYm,distribution{k});
                    [h(k,1),p(k,1)] = chi2gof(stYm,'CDF',pd);
                catch
                    h(k,1)=-Inf;
                    p(k,1)=-Inf;
                end
            end
    
            if find(h==0) == 1
                bestfit = 1;
            else
                [~,bestfit] = max(p);
            end
    
            try
                bestpdY = fitdist(stYm,distribution{bestfit});
            catch
                continue
            end
    
            for n = 1:ntf
    
                X = tfexpression(:,n);
                meanX = mean(X);
                maxX = max(X);
                minX = min(X);
                stX = (X - meanX)/(maxX - minX);
    
                Xmp = logical(stX >= mean(stX));
                Xmn = logical(stX <= mean(stX));
    
                stYp = stYm(Xmp);
                stYn = stYm(Xmn);
    
                if size(stYp)<=size(stYn)
                    stYt=stYp;
                else
                    stYt=stYn;
                end
    
                if abs(mean(stYp) - mean(stYn)) >= 0.05
                    [h0,p0] = chi2gof(stYt,'CDF',bestpdY);
                    if  h0 == 1
                        pv(n,j) = p0;
                    else
                        pv(n,j) = 0;
                    end
                end
            end
    
        end
    end
    
    pv = [pv{:}];

    Net_st = cell(nnz(pv),3);
    
    [tfs_idx,gene_idx,pvs] = find(pv);
    pvs = -log(pvs);
    
    [Net_st{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_st{:,2}] = gene_names{gene_idx};
    [Net_st(:,3)] = num2cell(pvs);
    
    sNet_st = sortrows(Net_st,3,'descend');
    
end


