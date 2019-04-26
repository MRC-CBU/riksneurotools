
% Hacky function for thresholding connectivity matrices across groups
% Rik.Henson@mrc-cbu.cam.ac.uk   Nov 2018

function connectivity_stats(S);

try
    CM = S.CM;
catch
    error('Need to pass a Nsubject x Nroi x Nroi connectivity matrix in S.CM')
end

try
    X = S.X;
catch
    error('Need to pass a design matrix in S.X')
end

try
    contrasts = S.contrasts;
catch
    error('Need to pass a contrast matrix in S.contrasts')
end

try
    Grp = S.Grp;
catch
    error('Need to pass a Nsubject matrix indicating which group is subject belongs to')
end

try 
    GrpPlot = S.GrpPlot; 
catch
    GrpPlot = [];
%    GrpPlot = eye(length(unique(Grp))); 
end

try
    GrpPlotName = S.GrpPlotName;
catch
    for g=1:size(GrpPlot,2)
        GrpPlotName{g} = sprintf('Grp%d',g);
    end
end

try
    rnam = S.rnam;
catch
    for r=1:size(CM,2)
        rnam{r} = sprintf('ROI%d',r);
    end
end


try CorType = S.CorType; catch CorType = 'Unknown', end

try FWE = S.FWE; catch FWE = 1, end
try FDR = S.FDR; catch FDR = 0, end
try Unc = S.Unc; catch Unc = 0, end
try Alpha = S.Alpha; catch Alpha = 0.05, end

try isDcor = S.isDcor; catch isDcor = 0, end
try MeanRegress = S.MeanRegress; catch MeanRegress = 0, end

try PlotTmap = S.PlotTmap; catch PlotTmap = 1, end
try PlotPmap = S.PlotPmap; catch PlotPmap = 1, end



%% Mean Regression

if MeanRegress
    CMmr = NaN(size(CM));
    meanCM = squeeze(mean(mean(CM,3),2));
    mX = [meanCM ones(size(meanCM))];
    XpX = mX*pinv(mX);
    for r = 1:size(CM,2);
        y = squeeze(CM(:,r,:));
        CMmr(:,r,:) = y - XpX*y;
        CMmr(:,r,:) = CMmr(:,r,:) + repmat(mean(CM(:,r,:),1),[size(CM,1) 1 1]); % Adding back mean according to Linda;
    end
    CM = CMmr;
    CorType = ['MeanReg' CorType];
end

% Linda's faster code but need to reshape results back into matrix form
% CM = CM;
% meanCM = squeeze(mean(mean(CM,3),2));
% beta=[meanCM ones(size(meanCM))]\squeeze(CM(:,:));
% CMmr(:,:)=squeeze(CM(:,:))-[meanCM ones(size(meanCM))]*beta;
% CMmr(:,:)=squeeze(CMmr(:,:))+repmat(squeeze(nanmean(CM(:,:),1)),[length(meanCM) 1]);

% Mean unchanged
%figure; imagesc(squeeze(mean(CMmr,1))); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('CMmr(PreW=%d,Mean=%d)',PreWhiten,MeanVox)); %caxis([-1 1]); 
%set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);


%% Plot averages and T-stats for combinations of Groups

for g = 1:size(GrpPlot,1)
    ind = find(ismember(Grp,find(GrpPlot(g,:))));
    tmp = CM(ind,:,:);    
    if ~isDcor        
        TransData = atanh(tmp);  % Fisher transform (leading diagonals of 1 now Inf, but removed below)
    else
        TransData = log(CM(ind,:,:)+1); % Correct for bounded Dcor 0-1? (0.001 to stop log(0)=-Inf)
    end
    
    figure; imagesc(squeeze(mean(tmp,1))); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Mean %s:%s',CorType,GrpPlotName{g})); caxis([0 1]);
    set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[])
    
    %% plot the T-stats across all subjects (vs 0)
    
    CMTmap = squeeze(mean(TransData,1))./(squeeze(std(TransData,0,1))/sqrt(length(ind))); CMTmap(find(eye(size(CMTmap)))) = NaN;    
    CMPmap = t2p(CMTmap,size(CM,1)-1,0); CMPmap(find(eye(size(CMPmap)))) = 1;
    
    if PlotTmap
        figure; imagesc(CMTmap); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Tmap %s:%s',CorType,GrpPlotName{g})); 
        set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
    end
    if PlotPmap
        figure; imagesc(log10(CMPmap)); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Log10Pmap %s:%s',CorType,GrpPlotName{g})); 
        set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
    end
end


%% plot the T-stats comparing groups

contrasts = [contrasts zeros(size(contrasts,1),size(X,2)-size(contrasts,2))];

Nroi = size(CM,2);
Tmap = NaN(size(contrasts,1),Nroi,Nroi);
Pmap = Tmap; Bmap = Tmap;
%Alpha = Alpha / size(contrasts,1)

% Could speed-up by using a multivariate GLM, but not too slow at moment
for c = 1:size(contrasts,1)
    for ri = 1:Nroi
        for rj = (ri+1):Nroi
            y = squeeze(CM(:,ri,rj));  
            if ~isDcor; y = atanh(y); else y = log(y+1); end
            [Tmap(c,ri,rj),F,Pmap(c,ri,rj),df,R2,cR2,B,r,aR2,iR2] = glm(y,X,contrasts(c,:)',-1);
            Bmap(c,ri,rj) = contrasts(c,:)*B;
            Bmap(c,rj,ri) = Bmap(c,ri,rj);
            %Tmap(c,rj,ri) = Tmap(c,ri,rj);
            %Pmap(c,rj,ri) = Pmap(c,ri,rj);
        end
    end

    figure; imagesc(squeeze(Bmap(c,:,:))); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Mean %s: %s',CorType,num2str(contrasts(c,:))));
    set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);

    if PlotTmap
        figure; imagesc(squeeze(Tmap(c,:,:))); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Tmap %s: %s',CorType,num2str(contrasts(c,:))));
        set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
    end
    
    uPmap = squeeze(Pmap(c,:,:));
    if PlotPmap
        figure; imagesc(log10(uPmap)); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Log10Pmap %s: %s',CorType,num2str(contrasts(c,:))));
        set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
    end
    
     
    if Unc
%        ThrPmap = uPmap;
%        ThrPmap(find(uPmap>=Alpha)) = 1;
        ThrPmap = squeeze(Tmap(c,:,:));
        ThrPmap(find(uPmap>=Alpha)) = NaN;
        
        if PlotPmap & ~isempty(find(uPmap<Alpha))
            figure; imagesc(ThrPmap); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('Unc %s: %s',CorType,num2str(contrasts(c,:))));
            set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);           
        end
        
        fprintf('\nContrast %s, Alpha = %8.7f, Unc\n',num2str(contrasts(c,:)),Alpha)
%        [ri,rj] = find(squeeze(ThrPmap)<1);
        [ri,rj] = find(squeeze(~isnan(ThrPmap)));
        for i = 1:length(ri)
            fprintf('%s (%d) --------- %s (%d) (T=%3.2f, p=%1.6f)\n',rnam{ri(i)},ri(i),rnam{rj(i)},rj(i),Tmap(c,ri(i),rj(i)),Pmap(c,ri(i),rj(i)))
        end
    end
    
    if FWE
        CorAlpha = Alpha / ((Nroi-1)*Nroi/2);
%        ThrPmap = uPmap;
%        ThrPmap(find(uPmap>=CorAlpha)) = 1;
        ThrPmap = squeeze(Tmap(c,:,:));
        ThrPmap(find(uPmap>=CorAlpha)) = NaN;
        
        if PlotPmap & ~isempty(find(uPmap<CorAlpha))
            figure; imagesc(ThrPmap); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('FWE %s: %s',CorType,num2str(contrasts(c,:))));
            set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);           
        end
        
        fprintf('\nContrast %s, Alpha = %8.7f, FWE\n',num2str(contrasts(c,:)),Alpha)
%        [ri,rj] = find(squeeze(ThrPmap)<1);
        [ri,rj] = find(squeeze(~isnan(ThrPmap)));
        for i = 1:length(ri)
            fprintf('%s (%d) --------- %s (%d) (T=%3.2f, p=%1.6f)\n',rnam{ri(i)},ri(i),rnam{rj(i)},rj(i),Tmap(c,ri(i),rj(i)),Pmap(c,ri(i),rj(i)))
        end
    end
    
    if FDR
%        ThrPmap = uPmap;       
        ThrPmap = squeeze(Tmap(c,:,:));
        ind = find(triu(ones(size(uPmap)),1));       
        [AboveThr, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(uPmap(ind),Alpha,'dep','no');
        AboveThr = find(AboveThr);

% Holm's method
%         [Ranked,Order] = sort(uPmap(ind),'ascend');
%         N = length(Ranked);
%         AboveThr = [];
%         for f = 1:length(Ranked)
%             if Ranked(f) >= Alpha/(N+1-f);
%                 break
%             else
%                 AboveThr(end+1) = Order(f);
%             end
%         end
        ThrPmap(setdiff([1:length(ThrPmap(:))],ind(AboveThr))) = NaN;
        
        if PlotPmap & ~isempty(AboveThr)
            figure; imagesc(ThrPmap); colormap(jet); set(gca,'FontSize',12); colorbar; title(sprintf('FDR %s: %s',CorType,num2str(contrasts(c,:))));
            set(gca,'dataAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
        end
        
        fprintf('\nContrast %s, Alpha = %8.7f, FDR\n',num2str(contrasts(c,:)),Alpha)
%        [ri,rj] = find(squeeze(ThrPmap)<1);
        [ri,rj] = find(squeeze(~isnan(ThrPmap)));       
        for i = 1:length(ri)
            fprintf('%s (%d) --------- %s (%d) (T=%3.2f, p=%1.6f)\n',rnam{ri(i)},ri(i),rnam{rj(i)},rj(i),Tmap(c,ri(i),rj(i)),Pmap(c,ri(i),rj(i)))
        end
    end
end


return


    
 