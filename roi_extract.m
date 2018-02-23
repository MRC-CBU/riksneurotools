function ROI = roi_extract(S);

% Extract ROI data from image volume(s), either 4D (e.g., resting state
% timecourses) or 3D (e.g., grey matter signal). Gets data from image files
% themselves; does not require SPM.mat. 
%
% ROIfiles  = cell array of cell arrays, ie set of files for each ROI mask
% Datafiles = cell array of cell arrays, ie set of files for each subject
% no_svd = 1 means do not do SVD (eg if data large), ie also calculate first eigenvector 
% output_raw = whether to save all volume x voxel raw data (as well as mean/median/SVD)
% zero_rel_tol = maximal proportion of voxels with no data (0 var or NaN) in order to include
% zero_abs_tol = minimum number of voxels needed in order to include
% uROIvals = cell array of vector of values that define ROI for each ROI mask
% mask_space = whether to use voxels from mask image (1), or data image (0, default)
%
% Original by Rik Henson Feb 2012.
% Extended by Jason Taylor, 2012.
% Mask_space option added by Linda Geerligs, May 2016
% Mask_space option fixed after ">" threshold added, Rik Henson,  July 2017

try ROIfiles = S.ROIfiles;
catch
    error('Must pass cell array of ROI mask files')
end

try Datafiles = S.Datafiles;
catch
    error('Must pass cell array of Data files')
end

try output_raw = S.output_raw;
catch
    output_raw = 0;
end

try no_svd = S.no_svd;
catch
    no_svd = 0;
end

try zero_rel_tol = S.zero_rel_tol;
catch
    zero_rel_tol = 0.5;
end

try zero_abs_tol = S.zero_abs_tol;
catch
    zero_abs_tol = 20;
end

try uROIvals = S.uROIvals;
catch
    uROIvals = {};
end

if isempty(uROIvals);
    for m=1:length(ROIfiles)
        ROIvals{m} = [];
    end
else
    if ~iscell(uROIvals)
        for m=1:length(ROIfiles)
            ROIvals{m} = uROIvals;
        end
    else
        for m=1:length(ROIfiles)
            try
                ROIvals{m} = uROIvals{m};
            catch
                ROIvals{m} = [];
            end
        end
    end
end

try mask_space = S.mask_space;
catch
    mask_space = 0;
    warning('Using voxels in data space (not mask space)')
end

do_overwrite = 0;
verbose      = 0;
svd_tol  = 1e-6;

ROI = struct(); 

for s=1:length(Datafiles)
    fprintf('Doing Subject %d/%d \n',s,length(Datafiles));
    
    fl = Datafiles{s}';
    VY = spm_vol(fl);
    VY = [VY{:}];
    
    [~,yXYZmm] = spm_read_vols(VY(1)); % Assume all data in same space!
    % Get inverse transform (assumes all data files are aligned)
    Yinv  = inv(VY(1).mat);
        
    nr=0;     
    for m=1:length(ROIfiles)
        fprintf('\tDoing Mask %d/%d \n',m,length(ROIfiles));

        VM = spm_vol(ROIfiles{m});
        Minv = inv(VM.mat);
        
        [YM,mXYZmm] = spm_read_vols(VM);
        ROIfstem = spm_str_manip(ROIfiles{m},'rt');
           
        % Transform ROI XYZ in mm to voxel indices in data:
        yXYZind = Yinv(1:3,1:3)*mXYZmm + repmat(Yinv(1:3,4),1,size(mXYZmm,2));
        % Transform data XYZ in mm to voxel indices in mask:
        mXYZind = Minv(1:3,1:3)*yXYZmm + repmat(Minv(1:3,4),1,size(yXYZmm,2));
        % Transform data XYZ in mm to voxel indices in data:
        yyXYZind = Yinv(1:3,1:3)*yXYZmm + repmat(Yinv(1:3,4),1,size(yXYZmm,2));
        
        YMdata = spm_get_data(VM,mXYZind); 

        if isempty(ROIvals{m})
            ROIvals{m} = setdiff(unique(YM(:)),0); 
        elseif ROIvals{m}(1) == '>'
            thr = str2num(ROIvals{m}(2:end));
            f = find(YM>thr);
            YM(f) = 1;
            f = find(YMdata>thr);
            YMdata(f) = 1;
            ROIvals{m} = [1];
        end
 
        fprintf('\t\tDoing ROI (/%d):',length(ROIvals{m}));  
        for r=1:length(ROIvals{m})
            fprintf('.%d',r);
            nr = nr+1;
                        
            ROI(nr,s).ROIfile  = ROIfstem;
            ROI(nr,s).svd_tol  = svd_tol;           
            ROI(nr,s).zero_rel_tol = zero_rel_tol;           
            
            if mask_space
                f = find(YM==ROIvals{m}(r));
                d = spm_get_data(VY,yXYZind(:,f)); 
                ROI(nr,s).XYZ       = mXYZmm(:,f);
                ROI(nr,s).XYZcentre = mean(mXYZmm(:,f),2);
            else
                f = find(YMdata == ROIvals{m}(r));
                d = spm_get_data(VY,yyXYZind(:,f));
                ROI(nr,s).XYZ       = yXYZmm(:,f);
                ROI(nr,s).XYZcentre = mean(yXYZmm(:,f),2);
            end
            
            Nvox = size(d,2);
            if verbose
                fprintf('Region %d (%s = %d): %d ',r,ROIfstem,ROIvals{m}(r),length(f));
            end
            ROI(nr,s).ROIval    = ROIvals{m}(r);
            ROI(nr,s).numvox    = Nvox;
            
            if output_raw
                ROI(nr,s).rawdata = d;
            end
            
            % Check for zero-variance voxels:
            zero_vox = (var(d)==0 | isnan(var(d)));
            zero_count = sum(zero_vox); 
            ROI(nr,s).nonzerovox = Nvox-zero_count;
            
            if (zero_count/Nvox > zero_rel_tol) | ((Nvox - zero_count) < zero_abs_tol)
                %if verbose
                    fprintf('(%d nonzero) voxels -- FAILED (%d percent)!\n',Nvox-zero_count,100*zero_rel_tol);
                %end
                ROI(nr,s).mean     = repmat(NaN,size(d,1),1);
                ROI(nr,s).median   = repmat(NaN,size(d,1),1);
                ROI(nr,s).svd      = repmat(NaN,size(d,1),1);
                ROI(nr,s).svd_vox  = repmat(NaN,size(d,2),1);
                ROI(nr,s).svd_pvar = NaN;
            else
                % Remove zero-variance voxels:
                f = setdiff(f,zero_vox);
               
                if mask_space
                    d = spm_get_data(VY,yXYZind(:,f)); 
                else
                    d = spm_get_data(VY,yyXYZind(:,f));
                end
                
                if verbose
                    fprintf('(%d nonzero) voxels\n',Nvox-zero_count);
                end
                
                % Remove voxel-wise mean if data from several volumes(?)
%                 if numel(d)>max(size(d))
%                     d = d-repmat(mean(d,1),size(d,1),1);
%                 end
                
                % MEAN/MEDIAN:
                ROI(nr,s).mean   = mean(d,2);
                ROI(nr,s).median = median(d,2);
                
                % SVD (only for timecourses) (note that singular vectors are L2 normalised by spm_en within spm_svd):
                if ~any(size(d)==1) & ~no_svd
                    [U,S,V] = spm_svd(d,svd_tol);
                    if isempty(S)
                        %if verbose
                            fprintf('..SVD FAILED!\n');
                        %end
                        ROI(nr,s).svd      = repmat(NaN,size(d,1),1);
                        ROI(nr,s).svd_vox  = repmat(NaN,size(d,2),1);
                        ROI(nr,s).svd_pvar = NaN;
                        ROI(nr,s).svd_tol  = svd_tol;
                    else
                        ROI(nr,s).svd      = full(U(:,1));
                        ROI(nr,s).svd_vox  = full(V(:,1));
                        ROI(nr,s).svd_pvar = S(1)/sum(full(S(:)));
                        ROI(nr,s).svd_tol  = svd_tol;
                    end
                    %if isempty(S)
                    %    svd_tol = 1e-9; % could increase tolerance?
                    %    [U,S,V] = spm_svd(dd,svd_tol);
                    %end
                else
                    ROI(nr,s).svd      = NaN;
                    ROI(nr,s).svd_vox  = NaN;
                    ROI(nr,s).svd_pvar = NaN;
                    ROI(nr,s).svd_tol  = NaN;
                end
            end
        end
        fprintf('\n')
    end
end
