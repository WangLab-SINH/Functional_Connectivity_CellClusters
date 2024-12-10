function genePls(ResponseVar,PredictVar,ResponseVarNames,PredictVarNames,cur_responsevar,dim,bootnum,outdir)
mkdir(outdir)
index = find(strcmp(ResponseVarNames,cur_responsevar));
y = ResponseVar(:,index);
Y = zscore(y);
X = zscore(PredictVar);

%% Select the principal component number
comp_idx = zeros(1, 1000);
 parfor col = 1:1000
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim, 'cv',4);
    [~, comp_idx(1,col)] = min(MSE(2,2:end));
 end
 
unique_values = unique(comp_idx);
counts = histcounts(comp_idx, [unique_values, unique_values(end)+1]);

[~,id1] = max(counts);
ncomp = unique_values(1,id1);

%%
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,ncomp);

%%
bootbeta = zeros(size(X,2),bootnum);
rng(1234);
parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,stats]=plsregress(Xr,Yr,ncomp); %perform PLS for resampled data
    bootbeta(:,i) = BETAr(2:end,1);%store (ordered) weights from this bootstrap run
end

%%
BETA = BETA(2:end,1);
betasd = std(bootbeta');
temp1 = BETA./betasd';
[Z1,ind1] = sort(temp1,'descend');
BETA = BETA(ind1);
PLS1 = PredictVarNames(ind1);


fid1 = fopen([outdir,'PLS_gene_beta.csv'],'w');
fprintf(fid1,'%s, %s, %s\n', 'GeneName','BETA','Zscore');
for i=1:length(PredictVarNames)
    fprintf(fid1,'%s, %f, %f\n', PLS1{i},BETA(i),Z1(i));
end
fclose(fid1);
end