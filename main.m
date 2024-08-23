%% Original image
load('original_image_128.mat')
figure;
imagesc(original_image)
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Original Image');
xlabel('Position [x]');
ylabel('Position [y]');

%% Load Sinogram
load('128_18.mat')
figure;
imagesc(sino)
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('CT Scan Sinogram');
xlabel('Angle [\theta]');
ylabel('Position [x]');

%% General variables
[numProjections, numAngles] = size(sino);

angleInterval = [30, 180];
angles = linspace(angleInterval(1), angleInterval(2), numAngles-2);

gridSize = 128;
l = -92:92;

%% A2
angleSample = 1;
projectionAngle = angles(angleSample);
projectionAngleRad = deg2rad(projectionAngle);

p = sino(:, angleSample);

% % v1 naïve approach
% b30_naive = zeros(gridSize, gridSize);
% for y = -63:64
%   for x = -63:64
%     l_sample = round(x * cos(projectionAngleRad) + y * sin(projectionAngleRad));
%     [numDatapoints, z] = size(sino(:, angleSample));
%     indexMidpoint = (numDatapoints - 1) / 2;
% 
%     lIndex = indexMidpoint + l_sample;
%     b30_naive(y + 64, x + 64) = sino(lIndex, angleSample);
%   end
% end
% 
% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30_naive);
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('2D Back Projection b(x,y)');
% subtitle('Using naïve pixel looping');
% xlabel('Position [x]');
% ylabel('Position [y]');

% v2
[x,y] = meshgrid(linspace(-63, 63, gridSize));
l_cartesian = x * cos(projectionAngleRad) + y * sin(projectionAngleRad);
b30 = interp1(l, p, l_cartesian, 'linear', 0);

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('2D Back Projection b(x,y)');
subtitle('Using meshgrid');
xlabel('Position [x]');
ylabel('Position [y]');

%% A3
image = zeros(gridSize, gridSize);
for idx = 1:numel(angles)
  projectionAngle = angles(idx);
  projectionAngleRad = deg2rad(projectionAngle);
  
  forward_projection = sino(:, angleSample);
  
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  l_cartesian = x * cos(projectionAngleRad) + y * sin(projectionAngleRad);
  backProjection = interp1(l, forward_projection, l_cartesian, 'linear', 0);
  image = image + backProjection;
end

image = (pi / (2 * numAngles)) * image;

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Back Projection image f_b(x,y)');
xlabel('Position [x]');s
ylabel('Position [y]');

%% B2

image = zeros(gridSize, gridSize);

omega = linspace(-1, 1, numProjections);
rampFilter = abs(omega);

for idx = 1:numel(angles) 
    p = sino(:, idx);
    
    P = fftshift(fft(p));
    P_filtered = P .* rampFilter';
    p_filtered = ifft(ifftshift(P_filtered));

    projectionAngle = angles(idx);
    projectionAngleRad = deg2rad(projectionAngle);

    [X, Y] = meshgrid(linspace(-64, 63, gridSize), linspace(-64, 63, gridSize));
    s = X * cos(projectionAngleRad) - Y * sin(projectionAngleRad);

    backprojected = interp1(l, p_filtered, s, 'linear', 0);
    
    image = image + backprojected;
end

% Normalize the reconstructed image
image = (pi / (2 * numAngles)) * image;

% Display the reconstructed image
figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Filtered Back Projection Image');
subtitle('Using ramp filter \omega \in [-1, 1]')
xlabel('Position [x]');
ylabel('Position [y]');

%% B3