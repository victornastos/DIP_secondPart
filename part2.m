% Load image
I = imread('im2.jpg');
% Convert to grayscale
I = rgb2gray(I);
%I = I(1:8:end,1:8:end);
imshow(I,'InitialMagnification','fit');
%resize of image
N = 10;
[rows, columns, channel] = size(I);
numberOfRows = round(rows/N);
numberOfColumns = round(columns/N);
image = imresize(I, [numberOfRows numberOfColumns]);

myCorners = myDetectHarrisFeatures(image);

function  corners = myDetectHarrisFeatures(Im)
%k = 0.04;
threshold = 0.01;
sigma = 1;
halfwid = sigma * 3;

[xx, yy] = meshgrid(-halfwid:halfwid, -halfwid:halfwid);

Gxy = exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

Gx = xx .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));
Gy = yy .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

numOfRows = size(Im, 1);
numOfColumns = size(Im, 2);

% 1) Compute x and y derivatives of image
Ix = filter2(Gx, Im);
Iy = filter2(Gy, Im);
size(Ix)

% 2) Compute products of derivatives at every pixel
Ix2 = Ix .^ 2;
Iy2 = Iy .^ 2;
Ixy = Ix .* Iy;

% 3)Compute the sums of the products of derivatives at each pixel
Sx2 = filter2(Gxy, Ix2);
Sy2 = filter2(Gxy, Iy2);
Sxy = filter2(Gxy, Ixy);

im = zeros(numOfRows, numOfColumns);

Rmax=0;
for i=1:numOfRows
    for j=1:numOfColumns
        M=im2double([Sx2(i,j) Sxy(i,j);Sxy(i,j) Sy2(i,j)]);
        R(i,j)=det(M)-0.01*(trace(M))^2;
        if(R(i,j)>Rmax)
            Rmax=R(i,j);
        end
    end
end

nx=size(R);
for i = 2:nx(1)-1
    for j = 2:nx(2)-1
        if(R(i,j)>threshold*Rmax &&R(i,j)> R(i-1,j-1) && R(i,j) > R(i-1,j) && R(i,j) > R(i-1,j+1) && R(i,j) > R(i,j-1) && R(i,j) > R(i,j+1) && R(i,j) > R(i+1,j-1) && R(i,j) > R(i+1,j) &&(( R(i,j) > R(i+1,j+1))|| R(i,j)<0))
            result(i,j) = 1;
        end
    end
end
[posc, posr] = find(result==1);
corners = [posc, posr];

imshow(Im,'InitialMagnification','fit');
hold on;
plot(posr,posc,'rs','MarkerSize',5);
end