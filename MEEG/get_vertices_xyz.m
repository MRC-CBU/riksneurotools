function Ip = get_vertices_xyz(D,xyz,rad,Nm,Pflag,remove_vert_thr,mods);

if nargin < 6
    remove_vert_thr = 0;
elseif nargin == 6
    error('Need to pass modalities if pass remove_vert_thr')
end

Np = size(xyz,2);

vert  = D.inv{D.val}.mesh.tess_mni.vert;

%% Get source vertices
Ip = cell(1,Np);
for p  = 1:Np
    Dp = sum([vert(:,1) - xyz(1,p), ...
        vert(:,2) - xyz(2,p), ...
        vert(:,3) - xyz(3,p)].^2,2);
    
    % nearest mesh points
    %--------------------------------------------------------------
    Ip{p} = find(Dp < rad^2)';
    
    if length(Ip{p}) < Nm, error('Insufficient vertices'); end
end

% Remove vertices in more than one ROI (not any?)
allIp = cat(2,Ip{:});
for p = 1:Np
    Ip{p} = setdiff(Ip{p},cat(2,Ip{setdiff(1:Np,p)}));
    if length(Ip{p}) < Nm, error('Insufficient vertices'); end
end


% Keep vertices with lowest gain correlation with vertices in other ROIs
if remove_vert_thr > 0
    
    forw = load(fullfile(D.path,D.inv{D.val}.gainmat)); % get label order in L
            
    % Here chose vertices with least overlap across all modalities - could choose different vertices for different modalities, but then Ip would have to be function of "m" too
    Lp = cell(1,Np); 
    for p = 1:Np
        Lp{p} = [];
        for m = 1:length(mods)
            Ic{m}   = indchantype(D, mods{m}, 'GOOD');
            IcL     = spm_match_str(forw.label, D.chanlabels(Ic{m})); % match order in L            
            Lp{p}   = [Lp{p}; forw.G(IcL,Ip{p})];
        end
    end
    
    for p = 1:Np
        oLp = cat(2,Lp{setdiff([1:Np],p)});
        c = [];
        for v = 1:size(Lp{p},2)
            c(v,:) = corr(Lp{p}(:,v),oLp);
        end
        cc = mean(abs(c),2);
        [~,ci] = sort(cc,'ascend');
        nc = round(length(cc)*remove_vert_thr);
        
        Ip{p} = Ip{p}(ci(1:nc));
    end
end

if Pflag
    figure
    Z = zeros(1,length(vert));
    for p=1:Np
        Col = Z; Col(Ip{p})=1;
        subplot(ceil(sqrt(Np)),ceil(sqrt(Np)),p)
        h = patch('vertices',D.inv{D.val}.mesh.tess_mni.vert,'faces',D.inv{D.val}.mesh.tess_mni.face,'FaceVertexCdata',Col','FaceColor','flat');
        axis image off
        view(90,-90)
        rotate3d on
        title(mat2str(xyz(:,p)))
    end
end
