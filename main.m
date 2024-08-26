clear all
%% Load Sinogram


% load('128_18.mat') % Part A
load('128_60.mat') % C3


% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), sino)
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('CT Scan Sinogram');
% xlabel('Angle [\theta]');
% ylabel('Position [x]');

%% General variables
[numProjections, numAngles] = size(sino);

angleInterval = [10, 180];
angles = linspace(angleInterval(1), angleInterval(2), numAngles);
angles

gridSize = 128;
l = -92:92;

%% A2
angleSample = 3;
projectionAngle = 30;
projectionAngleRad = deg2rad(projectionAngle);
disp(['Using angle ' num2str(projectionAngle) ' in rad ' projectionAngleRad]);
disp(['Using angle index ' num2str(angleSample) ' corresponding to ' num2str(angles(angleSample))])

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

% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30);
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('2D Back Projection b(x,y)');
% subtitle('Using meshgrid');
% xlabel('Position [x]');
% ylabel('Position [y]');

%% A3
image_a3 = zeros(gridSize, gridSize);
for idx = 1:numel(angles)
  projectionAngle = angles(idx);
  projectionAngleRad = deg2rad(projectionAngle);
  
  forward_projection = sino(:, angleSample);
  
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  l_cartesian = x * cos(projectionAngleRad) + y * sin(projectionAngleRad);
  backProjection = interp1(l, forward_projection, l_cartesian, 'linear', 0);
  image_a3 = image_a3 + backProjection;
end

image_a3 = (pi / (2 * numAngles)) * image_a3;

% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image_a3);
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('Back Projection image f_b(x,y)');
% xlabel('Position [x]');
% ylabel('Position [y]');

%% B2

image_b2 = zeros(gridSize, gridSize);

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
    
    image_b2 = image_b2 + backprojected;
end

% Normalize the reconstructed image
image_b2 = (pi / (2 * numAngles)) * image_b2;

% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image_b2);
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('Filtered Back Projection Image');
% subtitle('Using ramp filter \omega \in [-1, 1]')
% xlabel('Position [x]');
% ylabel('Position [y]');

%% C1
% Original image
load('original_image_128.mat')
% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), original_image)
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('Original Image');
% xlabel('Position [x]');
% ylabel('Position [y]');

mse_a3 = immse(image_a3, original_image);
mse_b2 = immse(image_b2, original_image);

disp(['MSE for A3: ' num2str(mse_a3)]);
disp(['MSE for B2: ' num2str(mse_b2)]);

%% C2
controlImage = iradon(sino, angles);
imageSize = size(controlImage);

% slice away outer rows
controlImage = controlImage(2:imageSize(1) - 1,:);
% slice away outer cols
controlImage = controlImage(:, 2:imageSize(2) - 1);

mse_iradon = immse(controlImage, original_image);
disp(['MSE for inverse Radon: ' num2str(mse_iradon)]);

% figure;
% imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), controlImage)
% axis image; % Ensure the aspect ratio is correct
% colormap gray; % Use a grayscale colormap
% colorbar; % Show a colorbar
% title('Original Image');
% subtitle('Reconstruction using iradon');
% xlabel('Position [x]');
% ylabel('Position [y]');
% 
