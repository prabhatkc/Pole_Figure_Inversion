
function plot_pf_layers(pf_3d, pf_strings, color)
s = size(pf_3d);

if (s(3)== 3)
    f1 = figure; set(gcf, 'Position', [110, 100, 1200, 400]);
    suptitle(pf_strings.title);
    subplot(1, 3, 1), imagesc(pf_3d(:, :, 1)), title(pf_strings.subtitle{1}), axis image off, colormap(color), colorbar;
    subplot(1, 3, 2), imagesc(pf_3d(:, :, 2)), title(pf_strings.subtitle{2}), axis image off, colormap(color), colorbar;
    subplot(1, 3, 3), imagesc(pf_3d(:, :, 3)), title(pf_strings.subtitle{3}), axis image off, colormap(color), colorbar;
elseif (s(3)== 4)
    f1 = figure; set(gcf, 'Position', [110, 100, 1200, 400]);
    suptitle(pf_strings.title);
    subplot(1, 4, 1), imagesc(pf_3d(:, :, 1)), title(pf_strings.subtitle{1}), axis image off, colormap(color), colorbar;
    subplot(1, 4, 2), imagesc(pf_3d(:, :, 2)), title(pf_strings.subtitle{2}), axis image off, colormap(color), colorbar;
    subplot(1, 4, 3), imagesc(pf_3d(:, :, 3)), title(pf_strings.subtitle{3}), axis image off, colormap(color), colorbar;
    subplot(1, 4, 4), imagesc(pf_3d(:, :, 4)), title(pf_strings.subtitle{4}), axis image off, colormap(color), colorbar;
end
end
