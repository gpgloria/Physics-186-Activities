clear all; close all;
%% Activity 2.1

%% Sinusoid along the x-direction, Amplitude ϵ [0,1], frequency is 4 cycles/cm

N = 200;                         % number of points                   
x = linspace(-2,2,N);            % generates a row vector of N points between -2 and 2
y = x;
[X,~] = meshgrid(x,y);           % 2-D grid coordinates of x and y

A_sin = 0.5*sin(2*pi*4*X)+0.5;   % Equation of the sinusoid 

figure(1); imshow(A_sin);        % 2-D Display
figure(2); mesh(x,y,A_sin);      % 3-D Display

%%  Grating, frequency is 5 line pairs/cm. A line pair is a strip of black (0) and white (1)
N = 200;                         % number of points 
x = linspace(-2,2,N);            % generates a row vector of N points between -2 and 2
y = x;
[X,~] = meshgrid(x,y);           % 2-D grid coordinates of x and y

R = 0.5*square(2*pi*5*X)+0.5;            % Grating using the square Function
A_grat = R;

figure(3); imshow(A_grat);       % 2-D Display
figure(4); mesh(x,y,A_grat);     % 3-D Display

%% Square Aperture, 1 cm x 1cm
N = 200;                   % number of points 
x = linspace(-2,2,N);      % generates a row vector of N points between -2 and 2
y = x;
[X,Y] = meshgrid(x,y);     % 2-D grid coordinates of x and y
 
A = zeros(size(X));        % 200x200 Matrix of element 0
B = zeros(size(X));        % 200x200 Matrix of element 0

% Finding the area of the square and equate to 1
A(((abs(X))<1)) = 1;       
B(((abs(Y))<1)) = 1;
A_sq = A+B;               % Combine

figure(5); imshow(A);   % 2-D Display
figure(6); mesh(x,y,A_sq); % 3-D Display

%%  Annulus, R(outer)= 2 cm, thickness of the ring is 0.25cm
N = 200;                          % number of points 
x = linspace(-2,2,N);             % generates a row vector of N points between -2 and 2
y = x;
[X,Y] = meshgrid(x,y);            % 2-D grid coordinates of x and y

R = sqrt(X.^2+Y.^2);              % Equation of a circle 
A = le(R,1.75);                   % logical array making R<1.75 = 1
B = le(R, 2);                     % logical array making R<2 = 1
A_annulus = B-A;                  % Difference of the 2 logical array

figure(5); imshow(A_annulus);
figure(6); mesh(x,y,A_annulus);

%%  Circular aperture with graded transmittance
N = 200;                   % number of points 
x = linspace(-2,2,N);      % generates a row vector of N points 
y = x;
[X,Y] = meshgrid(x,y);     % 2-D grid coordinates of x and y
sigma = 1;                 % Standard deviation of the Gaussian profile

R = sqrt(X.^2+Y.^2);       % Equation of the circle
A = zeros(size(R));        % Initialize a zero-element matrix

% Formula of a Gaussian Function with bias
g = (1/(sigma*sqrt(2*pi))).*exp(-0.5*(R.^2)/(sigma^2))-0.0875; 
A((R<1.75)) = 1;           % Make circle of radius < 1.75 be equal to 1
A((R>1.75)) = 0;           % Greater than 1.75, equate to 0
A_circ = A.*g;             % Combine Gaussian and circle

figure(7); imshow(A_circ);
figure(8); mesh(x,y,A_circ);

%% Olympic Rings
clear all ; close all;

% Initializations
N = 1024;                  % number of points
x = linspace(-17,17,N);    % generates a row vector of N points between -17 and 17  
y = x;
[X,Y] = meshgrid(x,y);     % 2-D grid of x and y

% Initialize all 3 layers with white background using ones() 
Rd = ones(N,N);                      
Gn = Rd;
Blu = Rd;

% Draw colored circles
Rt = 3; RC = 5; deg = 30;R_k = 4.3;
xt = Rt*cosd(deg); yt = Rt*sind(deg);

% Black Circle | Black = 0
R_bl = sqrt((X.^2) + (Y+Rt).^2);        % Equation of the circle
Rd((R_bl>R_k) & (R_bl<RC)) = 0;         % Red = 0
Gn((R_bl>R_k) & (R_bl<RC)) = 0;         % Green = 0
Blu((R_bl>R_k) & (R_bl<RC)) = 0;        % Blue = 0

% Yellow Circle | Yellow = Red + Green
R_yl =sqrt((X+2*xt).^2 + (Y-2*yt).^2);  % Equation of the circle
Blu(((R_yl>R_k) & (R_yl<RC))) = 0;      % Blue = 0

% Green Circle | Blue & Red = 0
R_gn = sqrt((X-2*xt).^2 + (Y-2*yt).^2); % Equation of the circle
Rd (((R_gn>R_k) & (R_gn <RC))) = 0;     % Red = 0
Blu (((R_gn>R_k) & (R_gn <RC))) = 0;    % Blue = 0

% Red Circle | Blue & Green = 0
R_rd = sqrt((X-4*xt).^2 + (Y+2*yt).^2); % Equation of the circle
Gn(((R_rd>R_k) & (R_rd<RC))) = 0;       % Green = 0
Blu(((R_rd>R_k) & (R_rd<RC))) = 0;      % Blue = 0

% Blue Circle | Red & Green = 0
R_blu = sqrt((X+4*xt).^2 + (Y+2*yt).^2); % Equation of the circle
Rd(((R_blu>R_k) & (R_blu<RC))) = 0;      % Red = 0
Gn(((R_blu>R_k) & (R_blu<RC))) = 0;      % Blue = 0

% Combining three layers and Displaying the logo
I(:,:,1) = Rd;
I(:,:,2) = Gn;
I(:,:,3) = Blu;

% Displaying the olympic logo
figure; image(I); axis equal;

%% Saving
%Let I be the image matrix
I = uint8(A_circ*255);
imwrite(I, "circular_graded_transmittance.jpg");
imwrite(I, "circular_graded_transmittance.bmp");
imwrite(I, "circular_graded_transmittance.png");
imwrite(I, "circular_graded_transmittance.tif");
