
% TO DO
% ENSURE EEG HAS AVERAGE REFERENCE!!!!!
% could average over trials after (only if induced effects)
% add spatial/temporal SVD of channel data? (not worth it)

clear

%subnum = [2 3 5 14 8 9 10 11 12 16 15 17 18 23 24 25] % Mapping from openneuro to Scott's numbers (original numbers)
subnum = [2 3 5 14 8 9 10 12 16 15 17 18 23 24 25] % Excluding Scott's number 11 whose debriefing was too poor
%subnum = [1 2 3 5 6 8 9 10 11 12 14 15 16 17 18 19 23 24 25] % 19 that Scott has for fMRI (not all used in paper though?)

ROIname = {'lEVC','rEVC','lOFA','rOFA','lFFA','rFFA'};

fmri_rois = 0; % whether use fMRI or spheres
   
if fmri_rois    % Option 1 - fMRI clusters
    Np  = length(ROIname);
    owd = '/imaging/henson/users/sl01/FacesData';
else            % Option 2 - Peaks+Sphere
    xyz = [
        -12    15   -38    36   -42    42
        -87   -90   -86   -86   -56   -52
        -6     3   -14   -10   -20   -16
        ];
    Rad = 16; % radius
    Np = size(xyz,2);
end

%----------------


mods = {'MEG','MEGPLANAR','EEG'};
Nmod = length(mods);

Nm   = 1; % Number of spatial modes of gain matrix per ROI

% Whether want sign (polarity) of gain vectors to match that of first ROI after SVD 
% (not sure whether makes much sense when Nm>1, but no harm)
% match_signs = zeros(1,Np); % just keep arbitrary sign of SVD
match_signs = ones(1,Np);  % force to match that of first ROI
% match_signs = [1 1 1 -1 -1 -1]; % example if had 3 left lateral sources and 3 right, where want opposite signs across hemispheres

ci   = [1:3]; %D.nconditions;  % Assume trial-averaged at moment
Ncon = length(ci);
Twin = [0 500]/1000; % in secs
conn = [0.5 0.5 -1]; % Contrast of conditions
Fwin = [150 190]/1000; % Timewindow for above contrast
%Fwin = [200 250]/1000;

Nsub = length(subnum);
Nmth = 4;  % Number of methods below

Han  = 1; % whether to apply Hanning

remove_vert_thr = 0.2; % If >0, then keeps that proportion of vertices with lowest overlap with other vertex ROIs

M      = cell(Nsub,Nmth);  % inverse operator
erp    = cell(Nsub,Nmth);
facscr = nan(Nsub,Np,Nmth);
R2     = nan(Nsub,Nmth); F = R2;
RL     = cell(Nmod,1); RC = RL;
for sub = 1:Nsub
    spmfile = sprintf('/imaging/henson/Wakeman/ds000117/derivatives/DCM_Run09/sub-%02d/meg/dcm_ready_maMceffdspmeeg_sub-%02d_ses-meg_task-facerecognition_run-01_proc-sss_meg',sub,sub)
    D = spm_eeg_load(spmfile);
    
    if fmri_rois
        for p = 1:Np
            MaskImages{p} = [fullfile(owd,sprintf('subject_%02d',subnum(sub)),'stat_concat_deb',['VOI_' ROIname{p} '_mask']) '.nii'];
        end
        
%        If want same ROI for all subjects (ie not subject specific p<.001 that Scott used):
%         owd     = '/imaging/henson/users/rh01/Collaborations/Scott/GroupStats/stat_concat_deb';
%         ROIname = {'lEVC_r-l','rEVC_l-r','lOFA_FWE','rOFA_FWE','lFFA_FWE','rFFA_FWE'};
%         for p = 1:Np
%             MaskImages{p} = [fullfile(owd,ROIname{p}) '.nii'];
%         end
        
        Ip = get_vertices_fmri(D,MaskImages,Nm,(sub==1),remove_vert_thr,mods); % Show clusters only if first subject
    else
        Ip = get_vertices_xyz(D,xyz,Rad,Nm,(sub==1),remove_vert_thr,mods);
    end

%% Get leadfields for each modality and source
    %[L,D] = spm_eeg_lgainmat(D);
    forw = load(fullfile(D.path,D.inv{D.val}.gainmat)); % get label order in L

    Ic = cell(Nmod,1);
    CL = cell(Nmod,1); % leadfields (after SVD)
    SL = CL;           % scaled leadfields
    Lscale = nan(1,Nmod);
    Nchn = nan(Nmod,1);
    for m = 1:Nmod

        Ic{m}  = indchantype(D, mods{m}, 'GOOD');
        IcL    = spm_match_str(forw.label, D.chanlabels(Ic{m})); % match order in L
        
        Nchn(m)  = length(Ic{m});
        %R{m} = speye(Nchn(m),Nchn(m));
        Lp = forw.G(IcL,Ip{1});
        
      
        fprintf('\nProp. Var. explained for %s by Nm=%d\n',mods{m},Nm)
        [U,ss] = spm_svd(Lp',0); ss = full(diag(ss));
        ss     = sum(ss(1:Nm).^2)/sum(ss.^2); fprintf('\tROI%d: %3.2f',1,full(ss));
        U      = U(:,1:Nm);
        UL     = Lp*U;
        CL{m}  = UL;
        if sub == 1
            RL{m} = CL{m}; % First ROI is reference for sign of SVD below 
            RC{m} = D.chanlabels(Ic{m});
        else
            if match_signs(1) ~= 0
                com_chans = spm_match_str(RC{m},D.chanlabels(Ic{m}));
                for n = 1:Nm
                    CL{m}(:,n) = match_signs(1)*UL(:,n)*sign(UL(:,n)'*RL{m}(com_chans,n));
                end
            end
        end

        for p = 2:Np
            Lp = forw.G(IcL,Ip{p});
            [U,ss]  = spm_svd(Lp',0); ss = full(diag(ss));
            ss      = sum(ss(1:Nm).^2)/sum(ss.^2); fprintf('\tROI%d: %3.2f',p,full(ss));
            U       = U(:,1:Nm);
            UL      = Lp*U;
            if match_signs(p) ~= 0
                com_chans = spm_match_str(RC{m},D.chanlabels(Ic{m}));
                for n = 1:Nm                   
                    UL(:,n) = match_signs(p)*UL(:,n)*sign(UL(:,n)'*RL{m}(com_chans,n));
                end
            end
            CL{m}   = [CL{m} UL];
        end
        
        % Scale leadfields to same L2 norm
        Lscale(m) = sqrt(trace(CL{m}*CL{m}')/Nchn(m));
        SL{m} = CL{m}/Lscale(m);
        %    figure,imagesc(corrcoef(full(SL{m}))),caxis([-1 1]),colorbar,title(mods{m}),set(gca,'XTickLabel',ROIname,'YTickLabel',ROIname)
    end
    
    CLscale = [];
    for m = 1:Nmod
        CLscale = [CLscale; ones(Nchn(m),Np)*Lscale(m)];
    end
    fprintf('\n')

%% get data   
    ss   = [indsample(D,Twin(1)):indsample(D,Twin(2))];
    fs   = [indsample(D,Fwin(1)):indsample(D,Fwin(2))]; [~,fs] = intersect(ss,fs);
    Nsam = length(ss);

    W  = sparse(1:Nsam,1:Nsam,spm_hanning(Nsam));
     
    y = []; sy = []; Qe = cell(1,Nmod); yscale = nan(1,Nmod);
    off = 0;
    for m = 1:Nmod

        % Baseline correct (since data lowpass filtered?)
        %bc = mean(d(:)); % Don't want to mean-correct each channel/condition separately?
        bc = nan(length(Ic{m}),Nsam,Ncon);
        for c = 1:Ncon
            bc(:,:,c) = repmat(mean(D(Ic{m},[1:find(D.time==0)],c),2),1,Nsam);
        end
        
        d = D(Ic{m},ss,ci);
 
        d = d - bc;  
        
        % Hanning
        if Han
            for c = 1:size(d,1)
                d(c,:,:) = W*squeeze(d(c,:,:));
            end
        end
        
        d = reshape(d,[Nchn(m) Nsam*Ncon]); % concatenate trials across columns
        y = [y; d];

        yscale(m) = sqrt(trace(d*d')/Nchn(m));  % could scale across all trials?
        d = d/yscale(m);
        sy = [sy; d];
        
        ind = [1:Nchn(m)] + off; % hacky!
        off = off + Nchn(m);

        Qe{1,m} = zeros(sum(Nchn)); 
        for ij = 1:length(ind)
            Qe{1,m}(ind(ij),ind(ij)) = 1;
        end
    end
    
    Cyscale = [];
    for m = 1:Nmod
        Cyscale = [Cyscale; ones(Nchn(m),Nsam*Ncon)*yscale(m)];
    end
    CL = full(cat(1,CL{:}));
    SL = full(cat(1,SL{:}));
%    figure,imagesc(corrcoef(full(SL))),caxis([-1 1]),colorbar,title('Fused'),set(gca,'XTickLabel',ROIname,'YTickLabel',ROIname)
    
    %% Simple pinv option
    
    %% Unscaled
    M{sub,1} = pinv(CL);
    dn = M{sub,1}*y;
        
    SSR = sum(var(full(y - CL*dn),0,2));
    SST = sum(var(full(y),0,2));
    R2(sub,1) = 100*(SST - SSR)/SST;

    d = nan(Np,Nsam*Ncon);
    for p = 1:Np
        d(p,:) = sum(dn([1:Nm]+(p-1)*Nm,:),1);
    end

    d = reshape(d,[Np Nsam Ncon]);
    erp{sub,1} = d;
    facscr(sub,:,1) = squeeze(mean(d(:,fs,:),2))*conn';   
    
    %% Scaled
    M{sub,2} = pinv(SL);
    dn = M{sub,2}*sy;
    
    SSR = sum(var(full(sy - SL*dn),0,2));
    SST = sum(var(full(sy),0,2));
%     SSR = sum(var(full(sy.*Cyscale - (SL.*CLscale)*dn),0,2));
%     SST = sum(var(full(sy.*Cyscale),0,2));
    R2(sub,2) = 100*(SST - SSR)/SST;

    if Nm>1
        d = nan(Np,Nsam*Ncon);
        for p = 1:Np
            d(p,:) = sum(dn([1:Nm]+(p-1)*Nm,:),1);
        end
    else
        d = dn;
    end
     
    d = reshape(d,[Np Nsam Ncon]);
    erp{sub,2} = d;
    facscr(sub,:,2) = squeeze(mean(d(:,fs,:),2))*conn';
    
    
    %% REML ...
    yy = sy*sy';
    
    %% ... with one source component (and 3 sensor components)
     
    LQP = {}; LQPL = {};
    LQP{1}   = SL; LQPL{1}  = SL*SL';
    Q        = [Qe LQPL];
    
    %Q0          = exp(-2)*trace(yy)/size(sy,2);
    %[Cy,h,Ph,F] = spm_reml_sc(yy,[],Q,size(sy,2),-4,16,Q0); % too sparse!
    [Cy,h,Ph,F(sub,3)] = spm_reml(yy,[],Q,size(sy,2));
    
    %Cp    = sparse(0);
    LCp   = sparse(0);
    hp    = h(Nmod + (1:length(LQP)));
    for j = 1:length(LQP)
        %Cp  =  Cp + hp(j)*QP{j};
        LCp = LCp + hp(j)*LQP{j};
    end
    
    M{sub,3} = LCp'/Cy;    
    dn       = M{sub,3}*sy;
    
    SSR = sum(var(full(sy - SL*dn),0,2));
    SST = sum(var(full(sy),0,2));
    R2(sub,3) = 100*(SST - SSR)/SST;
    
    if Nm>1
        d = nan(Np,Nsam*Ncon);
        for p = 1:Np
            d(p,:) = sum(dn([1:Nm]+(p-1)*Nm,:),1);
        end
    else
        d = dn;
    end
    
    d = reshape(d,[Np Nsam Ncon]);
    erp{sub,3} = d;
    facscr(sub,:,3) = squeeze(mean(d(:,fs,:),2))*conn';

%     figure,
%     for p = 1:Np
%         subplot(Np/2,2,p)
%         plot(D.time(ss),squeeze(d(p,:,:)))
%         if p==1; title(sprintf('Sub %d, Mth %d: %s',sub,3,ROIname{p})); else, title(ROIname{p}); end
%         axis([D.time(ss(1)) D.time(ss(end)) min(d(:)) max(d(:))])
%     end

    %% ... with Np source components (and 3 sensor components) - virtually same as above when Nm=1
 
    LQP = {}; LQPL = {};
    for p = 1:Np
        QP{1,p} = zeros(Np*Nm);
        ind = [1:Nm] + (p-1)*Nm;
        QP{p}(ind,ind) = ones(Nm); % covariance too
        %QP{p}(ind,ind) = eye(Nm);
        LQP{1,p}  = SL*QP{p};
        LQPL{1,p} = LQP{p}*SL';
    end
    Q        = [Qe LQPL];
    
    [Cy,h,Ph,F(sub,4)] = spm_reml(yy,[],Q,size(sy,2));
    
    %Cp    = sparse(0);
    LCp   = sparse(0);
    hp    = h(Nmod + (1:length(LQP)));
    for j = 1:length(LQP)
        %Cp  =  Cp + hp(j)*QP{j};
        LCp = LCp + hp(j)*LQP{j};
    end
    
    M{sub,4} = LCp'/Cy;
    dn       = M{sub,4}*sy;
    
    SSR = sum(var(full(sy - SL*dn),0,2));
    SST = sum(var(full(sy),0,2));
    R2(sub,4) = 100*(SST - SSR)/SST;
    
    if Nm>1
        d = nan(Np,Nsam*Ncon);
        for p = 1:Np
            d(p,:) = sum(dn([1:Nm]+(p-1)*Nm,:),1);
        end
    else
        d = dn;
    end
    
    d = reshape(d,[Np Nsam Ncon]);
    erp{sub,4} = d;
    facscr(sub,:,4) = squeeze(mean(d(:,fs,:),2))*conn';

% conditional variance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%----------------------------------------------------------------------
%Cq    = Cp - sum(LCp.*M')';
    
end

R2
F

pval = [];
for mth = 1:Nmth
    allerps = [];
    for sub = 1:Nsub
        allerps(sub,:,:,:) = erp{sub,mth};
    end
    
    mallerps = squeeze(mean(allerps,1));
    
    figure,
    for p = 1:Np
        subplot(Np/2,2,p)
        plot(D.time(ss),squeeze(mallerps(p,:,:)))
        title(ROIname{p})
        axis([D.time(ss(1)) D.time(ss(end)) min(mallerps(:)) max(mallerps(:))])
    end
    subplot(Np/2,2,1), title(sprintf('Method %d: %s',mth,ROIname{1}))
    
    [~,~,pval(:,mth)] = onet(squeeze(facscr(:,:,mth)));
end
    
pval

return

%% Apply inverse of choice to create LFP

mth = 3;

%% !!NOTE not written for Nm>1 yet - for that, would have to call spm_eeg_montage twice, 
%% first with Np*Nm inverse operator, and then with a second montage to average Nm 
%% neighbouring modes to get just Np channels

for sub = 1:Nsub
    spmfile = sprintf('/imaging/henson/Wakeman/ds000117/derivatives/DCM_Run09/sub-%02d/meg/dcm_ready_maMceffdspmeeg_sub-%02d_ses-meg_task-facerecognition_run-01_proc-sss_meg',sub,sub)
    D = spm_eeg_load(spmfile);
    
    S = [];
    S.D = D;
    S.mode = 'write';
    S.montage.tra = M{sub,mth};
    S.montage.labelnew = ROIname;
    S.montage.labelorg = chanlabels(D,indchantype(D,mods));
    S.montage.chantypenew = repmat({'LFP'},1,Np);
    S.keepothers = 0;
    S.keepsensors = 0;
    D = spm_eeg_montage(S);
end


