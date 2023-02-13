function Ip = get_vertices_fmri(D,MaskImages,Nm,Pflag,remove_vert_thr,mods);

Np = length(MaskImages);

M  = D.inv{D.val}.mesh.tess_mni;

Smooth = 0;

Ip = cell(1,Np); allIp = [];
for p = 1:Np
    
    V = spm_vol(MaskImages{p});
    q = spm_read_vols(V);
    q = q > 0.5;  % !!depends on nature of mask
    
    %     q = spm_mesh_project(m.vert,struct('dat',double(q),'mat',V.mat),'nn');
    %     q = spm_mesh_smooth(struct('faces',double(m.face),'vertices',m.vert), q', Smooth);
    %     q = q .* (q > exp(-8));
    %     f = find(q)';
    
    I = ones(V.dim);
    v = V.mat\[M.vert, ones(size(M.vert,1),1)]'; v = v(1:3,:)';    
    [vor, dist] = spm_voronoi(I, v, 'd5711');
    f = unique(vor(find(q)))';
    
    Ip{p} = f;
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
    figure(Pflag)
    Z = zeros(1,size(M.vert,1));
    for p=1:Np
        Col = Z; Col(Ip{p})=1;
        subplot(ceil(sqrt(Np)),ceil(sqrt(Np)),p)
        h = patch('vertices',M.vert,'faces',M.face,'FaceVertexCdata',Col','FaceColor','flat');
        axis image off
        view(90,-90)
        rotate3d on
        title(spm_str_manip(MaskImages{p},'t'))
    end
end

return

N = nifti(MaskImages{1});
d = N.dat(:,:,:);
[x,y,z] = ind2sub(size(d),find(d));
z = mode(z);

d = size(N.dat);
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph)
h_ctx   = patch('vertices',m.vert,'faces',m.face,'EdgeColor','b','FaceColor','b');
hold on
f1 = N.dat(:,:,z);
M  = N.mat;
[x,y,z] = ndgrid(1:d(1),1:d(2),z);
x1 = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
y1 = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
z1 = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);

s  = surf(x1,y1,z1,f1);
set(s,'EdgeColor','none')
axis image off;
colormap('gray');
view(180,-80);
rotate3d on
drawnow
hold off