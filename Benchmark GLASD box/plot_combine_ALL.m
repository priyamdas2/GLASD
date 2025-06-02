% Read all images
img1 = imread('plot_ackley.png');
img2 = imread('plot_griewank.png');
img3 = imread('plot_rastrigin.png');
img4 = imread('plot_rosenbrock.png');
img5 = imread('plot_sumsquares.png');
img6 = imread('plot_legend_only.png');

% Create a very tight layout
figure;
t = tiledlayout(2, 3, ...
    'TileSpacing', 'none', ...
    'Padding', 'none');

% Display images
nexttile; imshow(img1); axis off;
nexttile; imshow(img2); axis off;
nexttile; imshow(img3); axis off;
nexttile; imshow(img4); axis off;
nexttile; imshow(img5); axis off;
nexttile; imshow(img6); axis off;

% Optional: Resize figure
set(gcf, 'Position', [100, 100, 1400, 750]);

% Save to file
exportgraphics(gcf, 'plot_combined_ALL.png', 'Resolution', 300);