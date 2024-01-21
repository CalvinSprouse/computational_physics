% script start
clear;

% load stm data
stm_mat = load("stm.txt");

% plot the stm mat as a heatmap
plt = imagesc(stm_mat);
set(gca, FontSize=14, FontName="Times New Roman")
cmap = colormap("parula");
clrbar = colorbar;

set(gca, XTickLabel=[])
set(gca, YTickLabel=[])
% set(clrbar, YTickLabel=[])

ylabel(clrbar, "\textbf{Height}", Interpreter="latex")
% ylabel("vertical position")
% xlabel("horizontal position")
title("\textbf{Scanning Tunnel Microscopy of unknown material}", Interpreter="latex")

exportgraphics(gcf, "stm_plt.png", Resolution=300, BackgroundColor="none")