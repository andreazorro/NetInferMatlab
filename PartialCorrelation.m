function [sNet_pcp, sNet_pcs] = PartialCorrelation(gene_names, regulators, expressiondata)

    warning ('off','all')

    tfs = find(ismember(gene_names,regulators));

    tfexpression = expressiondata(:,tfs);

    ntf = numel(tfs);
    ngenes = size(expressiondata,2);
    
    rho_pcp = zeros(ntf,ngenes);
    pv_pcp = zeros(ntf,ngenes);

    rho_pcs = zeros(ntf,ngenes);
    pv_pcs = zeros(ntf,ngenes);
    
    for i = 1:ngenes

        expredatap_2 = expressiondata(:,i);

        for j = 1:ntf   

            expredatap_1 = tfexpression(:,j); 
            tfexprem = tfexpression; 
            tfexprem(:,j) = [];
            if (tfs(ntf)~=i)
                tfexprem(:,logical(tfs == i)) = [];
            end

            [rho_pcp(j,i),pv_pcp(j,i)] = partialcorr(expredatap_1, expredatap_2,tfexprem,'Type','Pearson');
            [rho_pcs(j,i),pv_pcs(j,i)] = partialcorr(expredatap_1, expredatap_2,tfexprem,'Type','Spearman');

        end
    end

    rho_pcp(pv_pcp>0.05 | rho_pcp<0.5 | rho_pcp==1) = 0;
    rho_pcs(pv_pcs>0.05 | rho_pcs<0.5 | rho_pcs==1) = 0;

    rho_pcp(isnan(rho_pcp)) = 0;
    rho_pcs(isnan(rho_pcs)) = 0;
    
    rho_pcp = abs(rho_pcp);
    rho_pcs = abs(rho_pcs);
    
    Net_pcp = cell(nnz(rho_pcp),3);
    Net_pcs = cell(nnz(rho_pcs),3);

    [tfs_idx,gene_idx,rho] = find(rho_pcp);
       
    [Net_pcp{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_pcp{:,2}] = gene_names{gene_idx};
    [Net_pcp(:,3)] = num2cell(rho);

    clear tfs_idx gene_idx rho

    [tfs_idx,gene_idx,rho] = find(rho_pcs);
       
    [Net_pcs{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_pcs{:,2}] = gene_names{gene_idx};
    [Net_pcs(:,3)] = num2cell(rho);

    sNet_pcp = sortrows(Net_pcp,3,'descend');
    sNet_pcs = sortrows(Net_pcs,3,'descend');





