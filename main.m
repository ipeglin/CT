clear all
%% Load Sinogram

load('128_18.mat') % Part A
% load('128_60.mat') % C3
[numProjections, numAngles] = size(sino);

%% General variables
l = -92:92;
gridSize = 128;

% define projection angles
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
  % back projection naïve approach
  b30Naive = zeros(gridSize, gridSize);
  for y = -63:64
    for x = -63:64
      lCartesian = round(x * cos(projectionAngleRad) - y * sin(projectionAngleRad));
      [numDatapoints, z] = size(sino(:, angleSample));
      indexMidpoint = (numDatapoints - 1) / 2;

      lIndex = indexMidpoint + lCartesian;
      b30Naive(y + 64, x + 64) = sino(lIndex, angleSample);
    end
  end
  
  figure;
  imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), b30Naive);
  axis image; % Ensure the aspect ratio is correct
  colormap gray; % Use a grayscale colormap
  colorbar; % Show a colorbar
  title('2D Back Projection b(x,y)');
  subtitle('Using naïve pixel looping');
  xlabel('Position [x]');
  ylabel('Position [y]');
  figName = ['plots/bp_30deg_' num2str(numel(angles))];
  print('-deps', figName)
  
  % interpolation technique
  forwardProjections = sino(:, angleSample);
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  lCartesian = x * cos(projectionAngleRad) - y * sin(projectionAngleRad);
  b30 = interp1(l, forwardProjections, lCartesian, 'linear', 0);
  
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
imageA3 = zeros(gridSize, gridSize);
for idx = 1:numel(angles)
  projectionAngle = angles(idx);
  angleSample = find(angles == projectionAngle);
  projectionAngleRad = deg2rad(projectionAngle);
  
  % get back projection
  forwardProjection = sino(:, angleSample);
  
  % interpolate
  [x,y] = meshgrid(linspace(-63, 63, gridSize));
  lCartesian = x * cos(projectionAngleRad) - y * sin(projectionAngleRad);
  backProjection = interp1(l, forwardProjection, lCartesian, 'linear', 0);
  imageA3 = imageA3 + backProjection;
end

% normalise image
imageA3 = (pi / (2 * numAngles)) * imageA3;

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), imageA3);
axis image; % Ensure the aspect ratio is correct
colormap gray; % Use a grayscale colormap
colorbar; % Show a colorbar
title('Back Projection image f_b(x,y)');
xlabel('Position [x]');
ylabel('Position [y]');
figName = ['plots/bp_' num2str(numel(angles))];
print('-deps', figName)

%% B2

imageB2 = zeros(gridSize, gridSize);

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
    
    % forward projection in frequency domain
    P = fftshift(fft(forwardProjections));

    % apply filter
    pFiltered = P .* rampFilter';
    pFiltered = ifft(ifftshift(pFiltered));

    projectionAngle = angles(idx);
    projectionAngleRad = deg2rad(projectionAngle);

    [X, Y] = meshgrid(linspace(-64, 63, gridSize), linspace(-64, 63, gridSize));
    lCartesian = X * cos(projectionAngleRad) - Y * sin(projectionAngleRad);

    % interpolate and add
    backprojected = interp1(l, pFiltered, lCartesian, 'linear', 0);
    imageB2 = imageB2 + backprojected;
end

% normalize image
imageB2 = (pi / (2 * numAngles)) * imageB2;

figure;
imagesc(linspace(-63, 64, gridSize), linspace(-63, 64, gridSize), imageB2);
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

mseA3 = immse(imageA3, original_image);
mseB2 = immse(imageB2, original_image);

disp(['MSE for A3: ' num2str(mseA3)]);
disp(['MSE for B2: ' num2str(mseB2)]);

%% C2
controlImage = iradon(sino, angles);
imageSize = size(controlImage);

% slice away outer rows
controlImage = controlImage(2:imageSize(1) - 1,:);
% slice away outer cols
controlImage = controlImage(:, 2:imageSize(2) - 1);

mseIradon = immse(controlImage, original_image);
disp(['MSE for inverse Radon: ' num2str(mseIradon)]);

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