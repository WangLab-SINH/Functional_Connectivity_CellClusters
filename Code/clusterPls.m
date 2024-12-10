function clusterPls(ResponseVar,PredictVar,ResponseVarNames,PredictVarNames,cur_responsevar,dim,bootnum,outdir)
    index = find(strcmp(ResponseVarNames,cur_responsevar));
    y = ResponseVar(:,index);
    Y = zscore(y);
    X = zscore(PredictVar);
    
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    CUMPCTVAR=cumsum(PCTVAR(2,:));
    [R1,p1]=corr(XS,Y);
    disp(PCTVAR(2,:));
    disp(R1');
    %%
    %assess significance of PLS result
    pctvar_m2=zeros(bootnum,15);
    cumsum_matrix=zeros(bootnum,15);
    rng(1234);
    parfor j=1:bootnum
        order=randperm(size(Y,1));
        Yp=Y(order,:);
        [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
        pctvar_m2(j,:)=PCTVARr(2,:);
        cumsum_matrix(j,:)=cumsum(PCTVARr(2,:));
    end
    
    %%
    cum_pvalues = zeros(1, size(cumsum_matrix, 2));
    for col = 1:size(cumsum_matrix, 2)
        cum_pvalues(col) = sum(cumsum_matrix(:, col) > CUMPCTVAR(col))/bootnum;
    end
    
    pct_pvalues = zeros(1, size(pctvar_m2, 2));
    for col = 1:size(pctvar_m2, 2)
        pct_pvalues(col) = sum(pctvar_m2(:, col) > PCTVAR(2,col))/bootnum;
    end
    disp(pct_pvalues)
    if pct_pvalues(1,1) < 0.05
        mkdir(outdir)
        writematrix(cum_pvalues, [outdir,'PLS_comp_cum_pvalues.csv']);
        writematrix(pct_pvalues, [outdir,'PLS_comp_pct_pvalues.csv']);
        %% Plot the percent of variance explained in the response variable as a function of the number of components
        %
        figure('visible','off');
        plot(1:15,cumsum(100*PCTVAR(2,:)),'-bo');
        xlabel('Number of PLS components');
        ylabel('Percent Variance Explained in y');
        saveas(gcf, [outdir,'CUMPCTVAR.pdf']);
        
        %
        figure('visible','off');
        plot(1:15,100*PCTVAR(2,:));
        xlabel('Number of PLS components');
        ylabel('Percent Variance Explained in y');
        saveas(gcf, [outdir,'PCTVAR.pdf']);
        
        writematrix(100*PCTVAR(2,:), [outdir,'PCTVAR.csv']);
        writematrix(cumsum(100*PCTVAR(2,:)), [outdir,'CUMPCTVAR.csv']);
        %%
        geneindex=1:size(PredictVarNames,2);
        %% correcte gene weight
        comp_idx = 1;   % which component is significant
        
        [R1,p1]=corr(XS,Y);
        if R1(comp_idx, 1)<0
            stats.W(:,comp_idx)=-1*stats.W(:, comp_idx);
            XS(:, comp_idx)=-1*XS(:, comp_idx);
        end
        [PLS1w, x1] = sort(stats.W(:, comp_idx),'descend');
        PLS1ids=PredictVarNames(x1);
        geneindex1=geneindex(x1);
        PLS1_score=XS(:, comp_idx);
        
        PLS1weights = zeros(size(X,2),bootnum);
        
        % parp = parpool(20);
        rng(1234);
        parfor i=1:bootnum
            myresample = randsample(size(X,1),size(X,1),1);
            res(i,:)=myresample; %store resampling out of interest
            Xr=X(myresample,:); % define X for resampled regions
            Yr=Y(myresample,:); % define Y for resampled regions
            [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
            
            temp=stats.W(:, comp_idx);%extract PLS1 weights
            newW=temp(x1); %order the newly obtained weights the same way as initial PLS
            if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
                newW=-1*newW;
            end
            PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
        end
        PLS1sw = std(PLS1weights');
        temp1=PLS1w./PLS1sw';
        [Z1,ind1]=sort(temp1,'descend');
        PLS1=PLS1ids(ind1);
        geneindex1=geneindex1(ind1);
        fid1 = fopen([outdir,'PLS1_geneWeights_mean.csv'],'w');
        for i=1:length(PredictVarNames)
            fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
        end
        fclose(fid1);
        writematrix(PLS1_score, [outdir,'PLS1_score.csv']);
    else
        disp('The first principal component of PLS was not significant');
    end

end