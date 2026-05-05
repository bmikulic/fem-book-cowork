function book_style(fig)
% BOOK_STYLE  Apply consistent typography to a figure for textbook export.
%
%   BOOK_STYLE(fig) bumps every axis/tick/title in the given figure to
%   12 pt and every legend to 11 pt, sets axes LineWidth to 0.9, and
%   ensures a visible bounding box. Call immediately before
%   EXPORTGRAPHICS so the saved PDF/PNG inherits the styled state.
%
%   BOOK_STYLE() with no argument uses GCF.
%
%   Line widths inside the plot are deliberately left untouched; scripts
%   typically pass 'LineWidth' explicitly to PLOT, and those choices
%   should not be overwritten here.

if nargin < 1
    fig = gcf;
end

ax = findall(fig, 'Type', 'axes');
set(ax, 'FontSize', 12, 'LineWidth', 0.9, 'Box', 'on');
set(findall(fig, 'Type', 'text'),   'FontSize', 12);
set(findall(fig, 'Type', 'legend'), 'FontSize', 11);

end
