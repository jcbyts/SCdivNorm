function h=plotWeights(xy, wf, sz, cmin, cmax)
% h=plotWeights(xy, wf, sz, cmin, cmax)
    if nargin < 5
        cmax=nan;
    end
    
    if nargin < 4
        cmin = nan;
    end
    if nargin < 3 || any(isnan(sz))
        sz = 200;
    end
    
    clrs = getColors(wf,cmin,cmax);

    hold on;
    axis off;
    axis equal;
    set(gcf,'color','w');
    if numel(sz) ~= numel(wf)
        sz = repmat(sz, numel(wf),1);
    end
    
    for ii = 1:numel(wf)
        h(ii)=plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz(ii), ...
            'Color', clrs(ii,:), 'LineStyle', 'none');
    end
end
