function C = shearPlotCanvas(tabGroup, name, nRows, nCols, figPos, headerText)
%SHEARPLOTCANVAS One plot group (a P#/Q# "figure") as a window or as a tab.
%  C = shearPlotCanvas([], name, nRows, nCols, figPos) opens a classic figure holding an
%  nRows x nCols subplot grid - what plotShearSoftening/plotShearBonds have always done.
%  C = shearPlotCanvas(tg, name, ...) instead adds a tab called `name` to the uitabgroup `tg`
%  and fills it with an nRows x nCols uigridlayout of uiaxes, so the plots stay inside the
%  app window and no new window is opened. figPos is ignored in that case.
%  headerText (optional, '' for none) is the group title: sgtitle in a figure, a label row in a tab.
%
%  Fields: C.container (the figure or the tab), C.ax (nRows*nCols axes, filled row-major),
%          C.isFigure, C.finish(outDir, fileName) - saves a PNG when outDir is given and
%          the group is a real figure (a tab cannot be exported with exportgraphics).

    if nargin < 6, headerText = ''; end
    nAx = nRows * nCols;
    C.isFigure = isempty(tabGroup);

    if C.isFigure
        C.container = figure('Color', 'w', 'Position', figPos);
        C.ax = gobjects(nAx, 1);
        if nAx == 1
            C.ax(1) = axes('Parent', C.container);
        else
            for q = 1:nAx
                C.ax(q) = subplot(nRows, nCols, q, 'Parent', C.container);
            end
        end
        if ~isempty(headerText)
            sgtitle(C.container, headerText, 'FontSize', 12, 'FontWeight', 'bold');
        end
    else
        C.container = uitab(tabGroup, 'Title', name, 'BackgroundColor', 'w');
        if isempty(headerText)
            g = uigridlayout(C.container, [nRows, nCols]);
            g.RowHeight = repmat({'1x'}, 1, nRows);
            rowOffset = 0;
        else
            g = uigridlayout(C.container, [nRows + 1, nCols]);
            g.RowHeight = [{44}, repmat({'1x'}, 1, nRows)];
            rowOffset = 1;
            lbl = uilabel(g, 'Text', headerText, 'FontWeight', 'bold', 'WordWrap', 'on', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'center');
            lbl.Layout.Row = 1;
            lbl.Layout.Column = [1 nCols];
        end
        g.ColumnWidth = repmat({'1x'}, 1, nCols);
        g.Padding = [8 8 8 8];
        g.RowSpacing = 8;
        g.ColumnSpacing = 8;

        C.ax = gobjects(nAx, 1);
        for q = 1:nAx
            a = uiaxes(g);
            a.Layout.Row = rowOffset + ceil(q / nCols);
            a.Layout.Column = 1 + mod(q - 1, nCols);
            C.ax(q) = a;
        end
    end

    C.finish = @(outDir, fileName) finishCanvas(C, outDir, fileName);
end

function finishCanvas(C, outDir, fileName)
    drawnow;
    try
        set(findall(C.container, 'Type', 'axes'), 'Toolbar', []);
    catch
        % uiaxes toolbars cannot always be cleared; harmless either way.
    end
    if C.isFigure && ~isempty(outDir)
        exportgraphics(C.container, fullfile(outDir, [fileName '.png']), 'Resolution', 170);
    end
end
