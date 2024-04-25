function [sNet_pe,sNet_ke,sNet_sp] = Correlation(gene_names, regulators, expressiondata)

    warning ('off','all')
   
    tfs = find(ismember(gene_names,regulators));

    tfexpression = expressiondata(:,tfs);

    parfor i=1:size(tfexpression,2)

    [rho_pe(i,:),pv_pe(i,:)] = corr(tfexpression(:,i), expressiondata,'Type','Pearson');
    [rho_ke(i,:),pv_ke(i,:)] = corr(tfexpression(:,i), expressiondata,'Type','Kendall');
    [rho_sp(i,:),pv_sp(i,:)] = corr(tfexpression(:,i), expressiondata,'Type','Spearman');

    end

    rho_pe(pv_pe>0.05 | rho_pe<0.5 | rho_pe==1) = 0;
    rho_ke(pv_ke>0.05 | rho_ke<0.5 | rho_ke==1) = 0;
    rho_sp(pv_sp>0.05 | rho_sp<0.5 | rho_sp==1) = 0;

    rho_pe = abs(rho_pe);
    rho_ke = abs(rho_ke);
    rho_sp = abs(rho_sp);

    Net_pe = cell(nnz(rho_pe),3);
    Net_ke = cell(nnz(rho_ke),3);
    Net_sp = cell(nnz(rho_sp),3);

    [tfs_idx,gene_idx,rho] = find(rho_pe);
       
    [Net_pe{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_pe{:,2}] = gene_names{gene_idx};
    [Net_pe(:,3)] = num2cell(rho);

    clear tfs_idx gene_idx rho

    [tfs_idx,gene_idx,rho] = find(rho_ke);

    [Net_ke{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_ke{:,2}] = gene_names{gene_idx};
    [Net_ke(:,3)] = num2cell(rho);

    clear tfs_idx gene_idx rho

    [tfs_idx,gene_idx,rho] = find(rho_sp);

    [Net_sp{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_sp{:,2}] = gene_names{gene_idx};
    [Net_sp(:,3)] = num2cell(rho);

    sNet_pe = sortrows(Net_pe,3,'descend');
    sNet_ke = sortrows(Net_ke,3,'descend');
    sNet_sp = sortrows(Net_sp,3,'descend');
      
%end


