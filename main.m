clear all
%% Load Sinogram


load('128_18.mat') % Part A
% load('128_60.mat') % C3
[numProjections, numAngles] = size(sino);

%% General variables
l = -92:92;
gridSize = 128;

angleInterval = [10, 180];
angles = linspace(angleInterval(1), angleInterval(2), numAngles);

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), sino)
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('CT Scan Sinogram');
xlabel('Angle [\theta]');
ylabel('Position [x]');
figName = ['plots/sinogram_' num2str(numel(angles))];
print('-deps', figName)

%% A2
if numAngles == 18
  projectionAngle = 30;
  projectionAngleRad = deg2rad(projectionAngle);
  angleSample = find(angles == projectionAngle);
  disp(['Found: ' num2str(angleSample) ' in ']);
  
  disp(['Using angle ' num2str(projectionAngle) ' in rad ' projectionAngleRad]);
  disp(['Using angle index ' num2str(angleSample) ' corresponding to ' num2str(angles(angleSample))])
  
  % 
  % v1 naïve approach
  b30_naive = zeros(gridSize, gridSize);
  for y = -63:64
    for x = -63:64
      l_sample = round(x * cos(projectionAngleRad) - y * sin(projectionAngleRad));
      [numDatapoints, z] = size(sino(:, angleSample));
      indexMidpoint = (numDatapoints - 1) / 2;

      lIndex = indexMidpoint + l_sample;
      b30_naive(y + 64, x + 64) = sino(lIndex, angleSample);
    end
  end
  
  figure;
  imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30_naive);
  axis image; % Ensure the aspect ratio is correct
  colormap gray; % Use a grayscale colormap
  colorbar; % Show a colorbar
  title('2D Back Projection b(x,y)');
  subtitle('Using naïve pixel looping');
  xlabel('Position [x]');
  ylabel('Position [y]');
  figName = ['plots/bp_30deg_' num2str(numel(angles))];
  print('-deps', figName)
  
  % v2
  forwardProjections = sino(:, angleSample);
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  l_cartesian = x * cos(projectionAngleRad) - y * sin(projectionAngleRad);
  b30 = interp1(l, forwardProjections, l_cartesian, 'linear', 0);
  
  figure;
  imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30);
  axis image; % Ensure the aspect ratio is correct
  colormap gray; % Use a grayscale colormap
  colorbar; % Show a colorbar
  title('2D Back Projection b(x,y)');
  subtitle('Using meshgrid');
  xlabel('Position [x]');
  ylabel('Position [y]');
  figName = ['plots/bp_30deg_' num2str(numel(angles))];
  print('-deps', figName)
end

%% A3
image_a3 = zeros(gridSize, gridSize);
for idx = 1:numel(angles)
  projectionAngle = angles(idx);
  angleSample = find(angles == projectionAngle);
  projectionAngleRad = deg2rad(projectionAngle);
  
  forward_projection = sino(:, angleSample);
  
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  l_cartesian = x * cos(projectionAngleRad) - y * sin(projectionAngleRad);
  backProjection = interp1(l, forward_projection, l_cartesian, 'linear', 0);
  image_a3 = image_a3 + backProjection;
end

image_a3 = (pi / (2 * numAngles)) * image_a3;

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image_a3);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Back Projection image f_b(x,y)');
xlabel('Position [x]');
ylabel('Position [y]');
figName = ['plots/bp_' num2str(numel(angles))];
print('-deps', figName)

%% B2

image_b2 = zeros(gridSize, gridSize);

omega = linspace(-1, 1, numProjections);
rampFilter = abs(omega);

figure;
plot(omega, rampFilter);
title('Rampfilter for FBP');
xlabel('Angle Frequency [\omega]');
ylabel('Magnitude Response');
print -deps fbp_filter

for idx = 1:numel(angles) 
    forwardProjections = sino(:, idx);
    
    P = fftshift(fft(forwardProjections));
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

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), image_b2);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Filtered Back Projection Image');
subtitle('Using ramp filter \omega \in [-1, 1]')
xlabel('Position [x]');
ylabel('Position [y]');
figName = ['plots/fbp_' num2str(numel(angles))];
print('-deps', figName)

%% C1
% Original image
load('original_image_128.mat')
figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), original_image)
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Original Image');
xlabel('Position [x]');
ylabel('Position [y]');
print -deps originalImage

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

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), controlImage)
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Original Image');
subtitle('Reconstruction using iradon');
xlabel('Position [x]');
ylabel('Position [y]');
figName = ['plots/iradon_' num2str(numel(angles))];
print('-deps', figName)


%% D1-3
load('original_image_128.mat')

resolutions = [18 60];
nMax = 10;

results = zeros(nMax, numel(resolutions));  

for idx = 1:numel(resolutions)
  lambda = resolutions(idx);

  for n = 1:nMax
    imageFile = ['data/128_' num2str(lambda) '_N=' num2str(n) '.mat'];
    load(imageFile)
    mse = immse(result_image, original_image);
    disp(['MSSE (' imageFile '): ' num2str(mse)])
    results(n, idx) = mse;
  end
end

% generate LaTeX table content
disp(['MSE of ART for ' num2str(nMax) ' iterations...']);
[numSims, numLambdas] = size(results);
string = '';
for row = 1:numSims
  string = [string ' & $N=' num2str(row) '$'];
end
string = [string ' \\'];
disp(string)

for col = 1:numLambdas
  string = [ '$PN=' num2str(resolutions(col)) '$'];
  

  for row = 1:numSims
    string = [string ' & ' num2str(round(results(row, col), 4))];
  end
  string = [string ' \\'];
  disp(string);
end

%% Part E
load('original_image_128.mat')

lambdas = [0.1 1];
nMax = 10;

results = zeros(nMax, numel(lambdas));  

for idx = 1:numel(lambdas)
  lambda = lambdas(idx);
    
    for n = 1:nMax
      imageFile = ['data/128_18_lambda=' num2str(lambda) '_N=' num2str(n) '.mat'];
      load(imageFile)
      mse = immse(result_image, original_image);
      disp(['MSSE (' imageFile '): ' num2str(mse)])
      results(n, idx) = mse;
    end
    
end

% generate LaTeX table content
disp(['MSE of ART with 128_18 and lambda in [' num2str(lambdas) ']']);
[numSims, numLambdas] = size(results);
string = '';
disp('\toprule')
for row = 1:numSims
  string = [string ' & $N=' num2str(row) '$'];
end
string = [string ' \\'];
disp(string);
disp('\midrule');

for col = 1:numLambdas
  string = [ '$\lambda=' num2str(lambdas(col)) '$'];
  

  for row = 1:numSims
    string = [string ' & ' num2str(round(results(row, col), 4))];
  end
  string = [string ' \\'];
  disp(string);
end
disp('\bottomrule');