I = imread('im2.jpg');
%I = I(1:8:end,1:8:end);
%resize of image
N = 10;
[rows, columns, ch] = size(I);
numberOfRows = round(rows/N);
numberOfColumns = round(columns/N);
image = imresize(I, [numberOfRows numberOfColumns]);
thet = 54*(pi/180);

myrotateImage = myImgRotation(image, thet);

function rotImg = myImgRotation(Im, angle)
[m,n,p]=size(Im);

diagm = sqrt(m^2+n^2);
diagn = sqrt(m^2+n^2);

for a=1:diagm
   for b=1:diagn
      i = uint16((a-diagm/2)*cos(angle)+(b-diagn/2)*sin(angle)+m/2);
      j = uint16(-(a-diagm/2)*sin(angle)+(b-diagn/2)*cos(angle)+n/2);
      if i>0 && j>0 && i<=m && j<=n           
         Irotate(a,b,:)=Im(i,j,:);
      end
   end
end
rotImg = Irotate;
figure;
imshow(Irotate);
end