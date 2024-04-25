delete(gcp('nocreate'))
% parpool('Processes',8)

if exist('InferNets','dir')==0 
    mkdir InferNets
end

dirout = "InferNets/";

file_tf = dir('*.tfs');                                                     
regulators = importdata(file_tf.name);

files_exp = dir('*.tsv');                                                  


for i = 1:length(files_exp)

filename_exp = files_exp(i).name;

A = importdata(filename_exp);

expressiondata=A.data;
gene_names=A.textdata;

gene_names(:,~any(expressiondata)) = [];
expressiondata(:,~any(expressiondata)) = [];

ngenes = size(expressiondata,2);

t = zeros(1,5);

% ANOVA

tic

[sNet_a] = Anova(gene_names, regulators, expressiondata);

writecell(sNet_a,strcat(dirout,'Anova_',filename_exp),'FileType','text','Delimiter','tab')

t(1,1) = toc;

% Friedman

tic

[sNet_f] = Friedman(gene_names, regulators, expressiondata);

writecell(sNet_f,strcat(dirout,'Friedman_',filename_exp),'FileType','text','Delimiter','tab')

t(1,2) = toc;

% Statmodel

tic

[sNet_st] = Statmodel(gene_names, regulators, expressiondata);

writecell(sNet_st,strcat(dirout,'Statmodel_',filename_exp),'FileType','text','Delimiter','tab')

t(1,3) = toc;

% Partial Correlation

tic

[sNet_pcp, sNet_pcs] = PartialCorrelation(gene_names, regulators, expressiondata);

writecell(sNet_pcp,strcat(dirout,'PartialCorrelationPearson_',filename_exp),'FileType','text','Delimiter','tab')
writecell(sNet_pcs,strcat(dirout,'PartialCorrelationSpearman_',filename_exp),'FileType','text','Delimiter','tab')

t(1,4) = toc;

% Correlation

tic;

[sNet_pe,sNet_ke,sNet_sp] = Correlation(gene_names, regulators, expressiondata);

writecell(sNet_pe,strcat(dirout,'Pearson_',filename_exp),'FileType','text','Delimiter','tab')
writecell(sNet_ke,strcat(dirout,'Kendall_',filename_exp),'FileType','text','Delimiter','tab')
writecell(sNet_sp,strcat(dirout,'Spearman_',filename_exp),'FileType','text','Delimiter','tab')

t(1,5) = toc;

times = cell(2,6);
times(1,:) = {'NumberofGenes','Anova','Friedman','Statmodel','PartialCorrelation','Correlation'};
times(2,:) = num2cell([ngenes,t]);

writecell(times,strcat(dirout,'NGenesandTime_',filename_exp),'FileType','text','Delimiter','tab')
delete(gcp('nocreate'));

end

