function varargout = share_caxis(ax_list, shared_clim)
% share_caxis(ax_list)
% share_caxis(ax_list, shared_clim)
% clim = share_caxis(ax_list)
% clim = share_caxis(ax_list, shared_clim)

if ~exist('shared_clim', 'var')
    shared_clim = cell2mat(arrayfun(@(ax) caxis(ax), ax_list, 'uni', 0));
    shared_clim = [min(shared_clim(:,1)), max(shared_clim(:,2))];
end

arrayfun(@(ax) caxis(ax, shared_clim), ax_list);

if nargout == 1
    varargout{1} = shared_clim;
elseif nargout > 1
    error('Only 1 output allowed');
end
    
end