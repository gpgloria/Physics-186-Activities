%% Applied Physics 187: Activity 3
% Gelli Mae P. Gloria
% 2015-06435

clear; close all; clc;

%% 3.2 Histogram Manipulation by Backprojection
% Reading image and getting its grayscale using rgb2gray
img = imread('D:\Acaaaaads\App Physics 186\Trial4.jpg');
img = rgb2gray(img);
figure(1); imshow(img);

% PDF and CDF of the original image
figure(2);
PDF = histogram(img(:),(0:255));
PDF.Normalization = 'probability';
[dist, bin] = imhist(img);
CDF = cumsum(dist)/numel(img);
figure(3); plot(CDF); title("Cumulative distribution function")

% Getting the desired CDFS: Linear, Quadratic, Logarithmic, and Sigmoid
x = (0:255);
y = linspace(-2,2,256);
desiredCDF = {x/255 x.^2/(255^2) (log(x+1))/max((log(x+1))) (erf(y)+1)/max((erf(y)+1))};

% Plotting the desired CDFs and displaying modified images
for i = 1:length(desiredCDF)
newGS = interp1(desiredCDF{i}, x, CDF(img(:)+1));
Igraynew = reshape(newGS, size(img));
PDF1 = hist(Igraynew(:),(0:255))/numel(Igraynew);
CDF1 = cumsum(PDF1);

% PDF of the enhanced image
figure(3+i);
bar(PDF1); title("New PDF")

% CDF Plotting
figure(7+i);
plot(CDF1, 'k-', 'LineWidth', 2,'MarkerEdgeColor','k')
hold on 
plot(CDF,'b--', 'LineWidth', 1.2,'MarkerEdgeColor','k')
plot(desiredCDF{i}, 'r-.','LineWidth', 1.5,'MarkerEdgeColor','k')
hold off;
legend('New CDF', 'Original CDF', 'Desired CDF')

% Enhanced Image
figure(11+i);
imshow(uint8(Igraynew)); title("Modified Image")
end

%% 3.3 Contrast Sketching
close all; clc; 

% Read image and get its grayscale form
img = imread('D:\Acaaaaads\App Physics 186\Trial8.jpg');
img = rgb2gray(img);

% Get the CDF of the original image
[counts, ~] = imhist(img);
cdf = cumsum(counts);
cdf = cdf/numel(img);

% Contrast sketching using Equation 2
I_max = 75; %interp1(cdf, x, 0.6);
I_min = 13; % interp1(cdf, x, 0.05);
I_new = double(img - I_min)/(I_max - I_min);
PDF1 = hist(I_new(:),(0:255))/numel(I_new);

subplot(1,2,1); imshow(img); title("Original Image")
subplot(1,2,2); imshow(I_new);title('Stretched Gray Scale image');

%% 3.4 White Balance
%% CONTRAST STRETCHING
close all; clc; clear
% Read image
img = imread('D:\Acaaaaads\App Physics 186\Trial6.jpg');
% img = imrotate(img, 90);
x = (0:255);
figure(3); imshow(img);

% Three color channels of the image:
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
figure(1)
subplot(2,2,1); imshow(R);
subplot(2,2,2); imshow(G);
subplot(2,2,3); imshow(B);

% Get the CDF of the three channels
%RED:
[red_count, ~] = imhist(R);
cdf_r = cumsum(red_count);
cdf_r = cdf_r/numel(R);

%GREEN:
[gr_count, ~] = imhist(G);
cdf_g = cumsum(gr_count);
cdf_g = cdf_g/numel(G);

%BLUE:
[bl_count, ~] = imhist(B);
cdf_b = cumsum(bl_count);
cdf_b = cdf_b/numel(B);

% Contrast sketching using Equation 2
% RED
R_max = double(interp1(cdf_r,x, 0.8));
R_min = interp1(cdf_r,x, 0.05);
I_restored(:,:,1) = double(R - R_min)/(R_max - R_min);

% GREEN
G_max = double(interp1(cdf_g,x, 0.8));
G_min = interp1(cdf_g,x, 0.05);
I_restored(:,:,2) = double(G - G_min)/(G_max - G_min);

% BLUE
B_max = double(interp1(cdf_b,x, 0.8));
B_min = interp1(cdf_b,x, 0.05);
I_restored(:,:,3) = double(B - B_min)/(B_max - B_min);

figure(2); imshow(double(I_restored));title('White balancing using Contrast Stretching');
[dist, bin] = imhist(I_restored);
CDF = cumsum(dist)/numel(I_restored);
figure(4); plot(CDF); title("Cumulative distribution function")
%% Gray World Algorithm
% Get the average value of the three channels
Rave = (mean(R(:)));
Gave = (mean(G(:)));
Bave = (mean(B(:)));

% Gray World Algorithm
I_gw(:,:,1) = (R)/Rave;
I_gw(:,:,2) = (G)/Gave;
I_gw(:,:,3) = (B)/Bave;
figure(5); imshow(double(I_gw));
title('White balancing using Gray World Algorithm');

%% White Patch Algorithm
clear all; close all; clc;
img = imread('D:\Acaaaaads\App Physics 186\Trial6.jpg');
x = (0:255);
figure(3); imshow(img);
I = imcrop(img);
figure(1);imshow(I);

R_img = img(:,:, 1);
G_img = img(:,:, 2);
B_img = img(:,:, 3);

R = I(:,:, 1);
G = I(:,:, 2);
B = I(:,:, 3);
Rave = mean(R(:));
Gave = mean(G(:));
Bave = mean(B(:));

I_gw(:,:,1) = R_img/Rave;
I_gw(:,:,2) = G_img/Gave;
I_gw(:,:,3) = B_img/Bave;
figure(5); imshow(double(I_gw));title('White balancing using White Patch Algorithm');