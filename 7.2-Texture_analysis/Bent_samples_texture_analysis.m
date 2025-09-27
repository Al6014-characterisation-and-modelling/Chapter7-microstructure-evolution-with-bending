%% Cube Texture Analysis

%% Setup and Parameters
% Crystal symmetries
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};

% Plotting convention
setMTEXpref('xAxisDirection', 'east');
setMTEXpref('zAxisDirection', 'outOfPlane');

% Analysis parameters
cube_threshold = 0.3;
tol_angle = 25 * degree;

% File selection
fname1 = 'A_texture.ctf';
fname2 = 'B_texture_highres.ctf';
fname3 = 'C_texture Data.ctf';
fname4 = 'D_texture Data.ctf';
fname5 = 'G_sample Data.ctf';
fname6 = 'H_texture.ctf';

fname = fname6; % Change this to select different file

%% Data Import and Processing
% Load EBSD data
ebsd = EBSD.load(fname, CS, 'interface', 'ctf', 'convertSpatial2EulerReferenceFrame');

% Apply rotation
rot = rotation('Euler', 0*degree, 0*degree, 0*degree);
ebsd = rotate(ebsd, rot, 'keepXY'); 

% Define sample directions
TD = vector3d.Z;

% Get sample name
[~, sampleName, ~] = fileparts(fname);


%% OPTIONAL STEP - Marking a small area to analyze 
plot(ebsd);

%region = [200 -500 1100 500]; %[xmin ymin xmax-xmin, ymax-ymin] B sample
%region = [0 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] B sample small

%region = [50 -500 1100 500]; %[xmin ymin xmax-xmin, ymax-ymin] D sample
%region = [80 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] D sample small
%region = [390 -500 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] D sample small2

%region = [0 -270 690 270]; %[xmin ymin xmax-xmin, ymax-ymin] C sample small
%region = [0 -350 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] C sample small2


%region = [250 -750 1350 675]; %[xmin ymin xmax-xmin, ymax-ymin] G sample
%region = [450 -470 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] G sample small

%region = [750 -470 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] G sample small2

region = [0 -650 1350 675]; %[xmin ymin xmax-xmin, ymax-ymin] H sample
%region = [250 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] H sample small
%region = [480 -380 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] H sample small2


%region = [250 -760 1350 660]; %[xmin ymin xmax-xmin, ymax-ymin] G sample - 125microns

%rectangle('position',region,'edgecolor','r','linewidth',2)

rectangle('position',region,'edgecolor','k','linewidth',1);
cropped_ebsd = ebsd(inpolygon(ebsd,region));
%IPF_map(cropped_ebsd, 'Aluminium', TD)
IPF_map(cropped_ebsd, 'Aluminium', TD)

% Use cropped data for all subsequent analysis
ebsd = cropped_ebsd;

%% Generate IPF-Z Map (before any thresholding)
figure('Position', [50, 50, 900, 700]);
IPF_map(ebsd('Aluminium'), 'Aluminium', vector3d.Z);

% Remove automatic MTEX scale bar
h = findall(gca, 'Type', 'text');
for i = 1:length(h)
    if contains(h(i).String, 'μm') || contains(h(i).String, 'micron')
        delete(h(i));
    end
end
% Remove white rectangular scale bar
h_rect = findall(gca, 'Type', 'rectangle');
for i = 1:length(h_rect)
    if isequal(h_rect(i).FaceColor, [1 1 1]) || isequal(h_rect(i).FaceColor, 'white')
        delete(h_rect(i));
    end
end

print(gcf, fullfile(pwd, [sampleName, '_EBSD_Map.png']), '-dpng', '-r300');

%% Cube Texture Analysis

% Calculate grains
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree);

% Define Cube orientation and find cube pixels
CubicCS = ebsd('Aluminium').CS;
ori_cube = orientation.byMiller([0 0 1], [1 0 0], CubicCS);
ebsd_cube = ebsd('Aluminium').findByOrientation(ori_cube, tol_angle);

% Calculate cube fraction for each grain
grainVolumeFraction = arrayfun(@(g) sum(ismember(ebsd_cube.grainId, g.id)) / g.area, grains);

% Select cube grains
cube_grains = grains(grainVolumeFraction >= cube_threshold);

% Create cube grains EBSD data
cubePixels = ismember(ebsd.grainId, cube_grains.id);
ebsd_cubeGrains = ebsd(cubePixels);

%% Generate PNG Files

% Map 1: Cube Fraction Map
figure('Position', [100, 100, 800, 600]);
plot(grains, grainVolumeFraction);
axis tight; axis equal; axis off;

% Add some padding to ensure scale bar is visible
xlim_current = xlim;
ylim_current = ylim;
x_range = xlim_current(2) - xlim_current(1);
y_range = ylim_current(2) - ylim_current(1);
xlim([xlim_current(1) - 0.05*x_range, xlim_current(2) + 0.05*x_range]);
ylim([ylim_current(1) - 0.1*y_range, ylim_current(2) + 0.05*y_range]);

% White-to-red colormap
custom_colormap = [1, 1, 1; 1, 0.8, 0.8; 1, 0.6, 0.6; 1, 0.3, 0.3; 1, 0, 0];
colormap(interp1(linspace(0,1,5), custom_colormap, linspace(0,1,256)));

cb = colorbar;
cb.FontSize = 12;
cb.Label.String = 'Cube fraction per grain';
cb.Label.FontSize = 14;
cb.Label.FontWeight = 'bold';
caxis([0, 1]);

hold on;
plot(grains(grainVolumeFraction >= cube_threshold).boundary, 'linewidth', 2, 'color', 'blue');
plot(grains(grainVolumeFraction < cube_threshold).boundary, 'linewidth', 1, 'color', 'black');
hold off;

print(gcf, fullfile(pwd, [sampleName, '_Cube_Fraction_Map.png']), '-dpng', '-r300');

% Map 2: Cube Grains Only
if ~isempty(cube_grains) && ~isempty(ebsd_cubeGrains)
    figure('Position', [150, 150, 900, 700]);
    IPF_map(ebsd_cubeGrains, 'Aluminium', TD);
    
    % Remove automatic MTEX scale bar
    h = findall(gca, 'Type', 'text');
    for i = 1:length(h)
        if contains(h(i).String, 'μm') || contains(h(i).String, 'micron')
            delete(h(i));
        end
    end
    % Remove white rectangular scale bar
    h_rect = findall(gca, 'Type', 'rectangle');
    for i = 1:length(h_rect)
        if isequal(h_rect(i).FaceColor, [1 1 1]) || isequal(h_rect(i).FaceColor, 'white')
            delete(h_rect(i));
        end
    end
    
    print(gcf, fullfile(pwd, [sampleName, '_Cube_Grains_Only.png']), '-dpng', '-r300');
end

% Map 3: Cube Grains with Boundaries
if ~isempty(cube_grains) && ~isempty(ebsd_cubeGrains)
    figure('Position', [200, 200, 900, 700]);
    IPF_map(ebsd_cubeGrains, 'Aluminium', TD);
    
    % Remove automatic MTEX scale bar
    h = findall(gca, 'Type', 'text');
    for i = 1:length(h)
        if contains(h(i).String, 'μm') || contains(h(i).String, 'micron')
            delete(h(i));
        end
    end
    % Remove white rectangular scale bar
    h_rect = findall(gca, 'Type', 'rectangle');
    for i = 1:length(h_rect)
        if isequal(h_rect(i).FaceColor, [1 1 1]) || isequal(h_rect(i).FaceColor, 'white')
            delete(h_rect(i));
        end
    end
    
    hold on;
    plot(cube_grains.boundary, 'linewidth', 2, 'color', 'black');
    hold off;
    print(gcf, fullfile(pwd, [sampleName, '_Cube_Grains_with_Boundaries.png']), '-dpng', '-r300');
end

% Map 4: Statistics Table
fig_stats = figure('Position', [300, 300, 800, 600], 'Color', 'white');
axis off;

stats_text = {
    '================================================================================';
    '                           CUBE TEXTURE ANALYSIS RESULTS';
    '================================================================================';
    sprintf('Sample:                    %s', sampleName);
    '--------------------------------------------------------------------------------';
    'ANALYSIS PARAMETERS:';
    sprintf('  Cube threshold:          %.2f (%.0f%%)', cube_threshold, cube_threshold*100);
    sprintf('  Tolerance angle:         %.1f°', tol_angle/degree);
    '  Grain boundary angle:    10.0°';
    '--------------------------------------------------------------------------------';
    'GRAIN STATISTICS:';
    sprintf('  Total grains analyzed:   %d', length(grains));
    sprintf('  Cube grains found:       %d', length(cube_grains));
    sprintf('  Cube grain percentage:   %.1f%%', 100*length(cube_grains)/length(grains));
    sprintf('  Non-cube grains:         %d', length(grains) - length(cube_grains));
    '--------------------------------------------------------------------------------';
    'CUBE FRACTION DISTRIBUTION:';
    sprintf('  0.0-0.2 fraction:        %d grains', sum(grainVolumeFraction >= 0.0 & grainVolumeFraction < 0.2));
    sprintf('  0.2-0.4 fraction:        %d grains', sum(grainVolumeFraction >= 0.2 & grainVolumeFraction < 0.4));
    sprintf('  0.4-0.5 fraction:        %d grains', sum(grainVolumeFraction >= 0.4 & grainVolumeFraction < 0.5));
    sprintf('  0.5-0.7 fraction:        %d grains', sum(grainVolumeFraction >= 0.5 & grainVolumeFraction < 0.7));
    sprintf('  0.7-1.0 fraction:        %d grains', sum(grainVolumeFraction >= 0.7 & grainVolumeFraction <= 1.0));
    '================================================================================';
};

text(0.05, 0.95, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'Courier New', 'FontSize', 10, 'Color', 'black');

print(fig_stats, fullfile(pwd, [sampleName, '_Statistics_Table.png']), '-dpng', '-r300');
close(fig_stats);

%% Summary
fprintf('Analysis complete. Generated files:\n');
fprintf('  %s_EBSD_Map.png\n', sampleName);
fprintf('  %s_Cube_Fraction_Map.png\n', sampleName);
fprintf('  %s_Cube_Grains_Only.png\n', sampleName);
fprintf('  %s_Cube_Grains_with_Boundaries.png\n', sampleName);
fprintf('  %s_Statistics_Table.png\n', sampleName);


%%

%% Pixel deviation from pure Cube (0..1)
% 0  -> exact cube
% 1  -> >= tol_angle away from cube (clipped)

ebsd_al = ebsd('Aluminium');
if isempty(ebsd_al)
    warning('No Aluminium pixels found. Skipping cube deviation map.');
else
    % Ideal Cube for fcc Al: {001}<100>
    cs = ebsd_al.CS;
    ori_cube = orientation.byMiller([0 0 1],[1 0 0], cs);

    % Misorientation angle (degrees) of each pixel to Cube
    ang_deg = angle(ebsd_al.orientations, ori_cube) / degree;

    % Normalize to 0..1 with clipping at tol_angle
    tol_deg = tol_angle / degree;
    dev_al  = min(ang_deg, tol_deg) / tol_deg;   % 0..1

    % Plot: Cube deviation map (0..1)
    fig_dev = figure('Position',[120,120,900,700]);
    plot(ebsd_al, dev_al);
    axis equal tight off

    % Custom red→white colormap
    ncol = 256;
    cmap = [linspace(1,1,ncol)', linspace(0,1,ncol)', linspace(0,1,ncol)']; % [R G B]
    % This goes from red [1 0 0] at 0 to white [1 1 1] at 1
    colormap(cmap);

    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = 'Cube deviation (0 = red, 1 = white)';
    cb.Label.FontWeight = 'bold';


    % Optional: faint background for non-Al pixels (light gray)
    hold on
    nonAl = ebsd(~ismember(ebsd.mineral, {'Aluminium'}));
    if ~isempty(nonAl)
        plot(nonAl, 'FaceColor', [0.92 0.92 0.92], 'EdgeColor', 'none');
    end
    hold off

    % Remove MTEX auto scale bar without Unicode matching
    try
        % Newer MTEX tags colorbars; if not, this just does nothing
        delete(findall(gca,'Tag','mtexColorbar'));
    catch
    end

    print(fig_dev, fullfile(pwd, [sampleName, '_CubeDeviation_Map.png']), '-dpng','-r300');
end

%% Pixel deviation from pure Cube (0..1)
% 0  -> exact cube (red)
% 1  -> >= tol_angle away from cube (blue)

ebsd_al = ebsd('Aluminium');
if isempty(ebsd_al)
    warning('No Aluminium pixels found. Skipping cube deviation map.');
else
    % Ideal Cube for fcc Al: {001}<100>
    cs = ebsd_al.CS;
    ori_cube = orientation.byMiller([0 0 1],[1 0 0], cs);

    % Misorientation angle (degrees) of each pixel to Cube
    ang_deg = angle(ebsd_al.orientations, ori_cube) / degree;

    % Normalize to 0..1 with clipping at tol_angle
    tol_angle = 0.5
    tol_deg = tol_angle / degree;
    dev_al  = min(ang_deg, tol_deg) / tol_deg;   % 0..1

    % Plot: Cube deviation map (0..1)
    fig_dev = figure('Position',[120,120,900,700]);
    plot(ebsd_al, dev_al);
    axis equal tight off

    % Custom red→blue colormap
    ncol = 256;
    cmap = [linspace(1,0,ncol)', ... % R: 1 → 0
            zeros(ncol,1), ...       % G: 0 → 0
            linspace(0,1,ncol)'];    % B: 0 → 1
    colormap(cmap);

    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = 'Cube deviation (0 = red, 1 = blue)';
    cb.Label.FontWeight = 'bold';

    % Optional: faint background for non-Al pixels (light gray)
    hold on
    nonAl = ebsd(~ismember(ebsd.mineral, {'Aluminium'}));
    if ~isempty(nonAl)
        plot(nonAl, 'FaceColor', [0.92 0.92 0.92], 'EdgeColor', 'none');
    end
    hold off

    % Remove MTEX auto scale bar if present
    try
        delete(findall(gca,'Tag','mtexColorbar'));
    catch
    end

    % Save PNG
    print(fig_dev, fullfile(pwd, [sampleName, '_CubeDeviation_Map_bluered.png']), '-dpng','-r300');
end



%% Pixel-wise deviation from pure Cube (0..1)
% Compare tolerance = 10° vs 25°

ebsd_al = ebsd('Aluminium');
if isempty(ebsd_al)
    warning('No Aluminium pixels found. Skipping cube deviation map.');
else
    % Ideal Cube for fcc Al: {001}<100>
    cs = ebsd_al.CS;
    ori_cube = orientation.byMiller([0 0 1],[1 0 0], cs);

    % Misorientation angle (degrees) of each pixel to Cube
    ang_deg = angle(ebsd_al.orientations, ori_cube) / degree;

    % Loop over two tolerance values
    for tol_deg = [40, 50]

        % Normalise to 0..1 with clipping at current tolerance
        dev_al = min(ang_deg, tol_deg) / tol_deg;

        % Plot: Cube deviation map (0..1)
        fig_dev = figure('Position',[120,120,900,700]);
        plot(ebsd_al, dev_al);
        axis equal tight off

        % Custom red→blue colormap
        ncol = 256;
        cmap = [linspace(1,0,ncol)', ... % R: 1 → 0
                zeros(ncol,1), ...       % G: 0 → 0
                linspace(0,1,ncol)'];    % B: 0 → 1
        colormap(cmap);

        caxis([0 1]);
        cb = colorbar;
        cb.Label.String = sprintf('Cube deviation (0 = red, 1 = blue), tol = %d°', tol_deg);
        cb.Label.FontWeight = 'bold';

        % Optional: faint background for non-Al pixels (light gray)
        hold on
        nonAl = ebsd(~ismember(ebsd.mineral, {'Aluminium'}));
        if ~isempty(nonAl)
            plot(nonAl, 'FaceColor', [0.92 0.92 0.92], 'EdgeColor', 'none');
        end
        hold off

        % Save PNG with tolerance in filename
        print(fig_dev, fullfile(pwd, ...
            sprintf('%s_CubeDeviation_%ddeg.png', sampleName, tol_deg)), ...
            '-dpng','-r300');
    end
end
