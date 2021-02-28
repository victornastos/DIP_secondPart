% Load image
I = imread('im2.jpg');
% Convert to grayscale
I = rgb2gray(I);
%I = I(1:8:end,1:8:end);
%resize of image
N = 5;
[rows, columns, channel] = size(I);
numberOfRows = round(rows/N);
numberOfColumns = round(columns/N);
image = imresize(I, [numberOfRows numberOfColumns]);
%Gaussian filter
sigma = 3.75;
I_smooth = imgaussfilt(image, sigma);

% Edge filter - use edge()
I_BW = edge(I_smooth, 'Sobel');
figure
imshow(I_BW)
[H,L,res] = myHoughTransform(I_BW, 1, pi/180, 15);

function [H,L,res] = myHoughTransform(img_binary , Drho , Dtheta, n)

I = img_binary;
[rows, cols] = size(I);
%Rads to degrees
thetaDeg = Dtheta * 180/pi ; 
theta_maximum = 90;
rho_maximum = (sqrt(rows^2 + cols^2));
thetaScale = -theta_maximum:thetaDeg:theta_maximum - 1;
rhoScale = -rho_maximum:Drho:rho_maximum;

H = zeros(length(rhoScale), length(thetaScale));
%Hough transform
for row = 1:rows
    for col = 1:cols
        if I(row, col) > 0
            x = col - 1;
            y = row - 1;
            for theta_i = 1 : length(thetaScale)
                tempR = x*cos(thetaScale(theta_i) *pi/180) + y*sin(thetaScale(theta_i) * pi/180);
                tempR = round((tempR + rho_maximum)/Drho)+1;
                H(tempR,theta_i) = H(tempR,theta_i) + 1;

            end
        end
    end
end
figure
%print Hough transform
imshow(H,[],'XData',thetaScale,'YData', rhoScale,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;

%Hough Peaks
rhos = zeros(n, 1); 
thetas = zeros(n, 1); 
H_new = H;
H_padded= padarray(H_new, [1, 1], 'replicate'); 
[rows, cols] = size(H_new); 

for i = 2:rows-1 % to account for padding
    for j = 2:cols-1
        if any(find((H_padded(i-1:i+1, j-1:j+1) > H_padded(i,j)))) > 0 % if any of the neighbors are greater than center pixel
            H_new(i-1,j-1) = 0; % non-maximal suppression
        end
    end
end

for i = 1:n
    maxIdx = max(H_new(:)); % highest score
    [rhoMaxIdx, thetaMaxIdx] = find(H_new==maxIdx);
    rhos(i) = rhoMaxIdx(1); % add - 1 term to account for padding an extra row and column
    thetas(i) = thetaMaxIdx(1);
    H_new(rhoMaxIdx(1), thetaMaxIdx(1)) = 0; % clear the highest scoring cell, then move on
end

%array of the Peaks
L=[rhos,thetas];
plot(thetaScale(L(:,2)),rhoScale(L(:,1)),'s','color','red');
title('Image with hough peaks');
% figure
% imshow(I)
% hold on;
size(L,1);
figure
I = imread('im2.jpg');
imshow(I)
hold on
%Hough Lines
for i = 1:size(L,1)
    rho_Temp = rhoScale(L(i,1));
    the_Temp = thetaScale(L(i,2));
    if the_Temp == 0
        lineX1 = rho_Temp;
        lineX2 = rho_Temp;
        lineY1 = 1;
        lineY2 = size(I,1);
    else
        lineX1 = 1;
        lineX2 = size(I,2);
        lineY1 = (rho_Temp - lineX1*cos(the_Temp*pi/180)) / sin(the_Temp*pi/180);
        lineY2 = (rho_Temp - lineX2*cos(the_Temp*pi/180)) / sin(the_Temp*pi/180);
    end
    %Res
    k(i) = sqrt((lineY2-lineY1)^2 + (lineX2-lineX1)^2);

    plot(5.*[lineX1,lineX2],5.*[lineY1,lineY2],'r','LineWidth',2);
    title('Image with hough lines');
end
res= (size(x,1)*size(x,2))- ceil(sum(k));
end

