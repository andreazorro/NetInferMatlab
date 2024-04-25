function [sNet_a] = Anova(gene_names, regulators, expressiondata)

    warning ('off','all')
   
    tfs = find(ismember(gene_names,regulators));

    ngenes = size(expressiondata,2);
    ntf = size(tfs,2);

    n2_a=zeros(ntf,ngenes);
  
    tfexpression = expressiondata(:,tfs);

    for i = 1:ngenes

        expredatap_2 = expressiondata(:,i)';

        parfor j = 1:ntf   
            
            if i ~= tfs(1,j)
                
                expredatap_1 = tfexpression(:,j)'; 

                [p_a,SS_a] = anova2([expredatap_1;expredatap_2],1,'off');
                
                if p_a(1,2) < 0.05 
                    np_a=SS_a{3,2}/SS_a{5,2};   %SSA/SStotal
                else
                    np_a=0;
                end
               
                expredatan_1 = -expredatap_1;
                expredatan_2 = expredatap_2;

                [p_a,SS_a] = anova2([expredatan_1;expredatan_2],1,'off');
         
                if p_a(1,2) < 0.05 
                    nn_a=SS_a{3,2}/SS_a{5,2};  %SSA/SStotal
                else
                    nn_a=0;
                end
               
                if np_a>nn_a
                    n2_a(j,i)=np_a;
                else
                    n2_a(j,i)=-nn_a;
                end
            end
        end
    end

    n2_a=abs(n2_a);
   
    Net_a = cell(nnz(n2_a),3);
        
    [tfs_idx,gene_idx,n2s] = find(n2_a);
    
    [Net_a{:,1}] = gene_names{tfs(tfs_idx)};
    [Net_a{:,2}] = gene_names{gene_idx};
    [Net_a(:,3)] = num2cell(n2s);

    sNet_a = sortrows(Net_a,3,'descend');
      
end


