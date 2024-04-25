function [sNet_f] = Friedman(gene_names, regulators, expressiondata)

    warning ('off','all')
   
    tfs = find(ismember(gene_names,regulators));

    ngenes = size(expressiondata,2);
    ntf = size(tfs,2);

    n2_f=zeros(ntf,ngenes); 

    tfexpression = expressiondata(:,tfs);

    for i = 1:ngenes

        expredatap_2 = expressiondata(:,i)';

        parfor j = 1:ntf   
            
            if i ~= tfs(1,j)
                
                expredatap_1 = tfexpression(:,j)'; 

                [p_f,SS_f] = friedman([expredatap_1;expredatap_2],1,'off');

                if p_f < 0.05 
                    np_f=SS_f{2,2}/SS_f{4,2};   %SSA/SStotal
                else 
                    np_f=0;
                end

                expredatan_1 = -expredatap_1;
                expredatan_2 = expredatap_2;

                [p_f,SS_f] = friedman([expredatan_1;expredatan_2],1,'off')

                if p_f < 0.05 
                    nn_f=SS_f{2,2}/SS_f{4,2};  %SSA/SStotal
                else 
                    nn_f=0;
                end

                if np_f>nn_f
                    n2_f(j,i)=np_f;
                else
                    n2_f(j,i)=-nn_f;
                end
            end
        end
    end

    n2_f=abs(n2_f);
   
    Net_f = cell(nnz(n2_f),3);
        
    [tfs_idx,gene_idx,n2s] = find(n2_f);
    
    [Net_f{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_f{:,2}] = gene_names{gene_idx};
    [Net_f(:,3)] = num2cell(n2s);

    sNet_f = sortrows(Net_f,3,'descend');
    
end


