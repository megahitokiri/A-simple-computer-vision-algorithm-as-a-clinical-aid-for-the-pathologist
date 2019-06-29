Matlab Code
rgb = imread('R:\No.3.jpg');
im=rgb;
R=rgb(:,:,1);
G=rgb(:,:,2);
B=rgb(:,:,3);
rgb2=rgb2gray(rgb);
[x,y,z]=size(rgb);
HematoxEosin=128;
purple=0;
pink=0;
other=0;
for i = 1:x 
    for j = 1:y
        R = rgb(i,j,1);
        G = rgb(i,j,2);
        B = rgb(i,j,3);
        if (R >= 70) && (R <= 150)
            if (G >= 60) && (G <= 100)
                if (B >= 145) && (B <= 255)
                HematoxEosin=0;
                purple=purple+1;
                end
            end
        end
        
        if (R >= 200) && (R <= 255)
            if (G >= 0) && (G <= 160)
                if (B >= 150) && (B <= 225)
                HematoxEosin=255;
                pink=pink+1;
                end
            end
        end
        if (R >= 0) && (R <= 69)
        HematoxEosin=128;
        other=other+1;
        end
        if (R >= 151) && (R <= 199)
        HematoxEosin=128;
        other=other+1;
        end
        rgb2(i,j)=HematoxEosin;%%trinary image
    end
end
figure();
imshow(rgb2);
[~, threshold] = edge(rgb2, 'sobel');
fudgeFactor = .60;
 
E = entropyfilt(rgb2);
Eim = mat2gray(E);
rgb2 = Eim;
figure();
imshow(rgb2);
 
%%otsu trheshold method.
threshold = graythresh(rgb2);
bw = im2bw(rgb2,threshold);
BW1 = edge(bw,'sobel');
ratio=(pink/purple)*100;
display(pink);
display(purple);
BWoutline = bwperim(BW1);
Segout = rgb;
Segout(BWoutline) = 100;%100 es rojo y azul, 255 solo rojo
rgbgreenout=Segout;
rgbgreenout(:,:,2)=0;
 
%showing inner border
 
gris=rgb2gray(im);
imR=double(im(:,:,1));
imG=double(im(:,:,2));
imB=double(im(:,:,3));
imR2=(imB-imG-imR);
masc=(imR2>20);
imR2=imR2.*masc;
imR2=medfilt2(imR2);
imR2=imR2/255;
imR3=imadjust(imR2,[],[],1.8);
threshold = graythresh(imR3);
imR4 =im2bw(imR3,threshold);
 
rgbBlue = im;
rgb2Blue = imR4;
[~, threshold] = edge(rgb2Blue, 'sobel');
fudgeFactor = .60;
 
E = entropyfilt(rgb2Blue);
Eim = mat2gray(E);
 
rgb2Blue = Eim;
threshold = graythresh(rgb2Blue);
bwBlue = im2bw(rgb2Blue,threshold);
BW1Blue = edge(bwBlue,'sobel');
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
 
BW2Blue = imdilate(BW1Blue, [se90 se0]);
BWdfillBlue = imfill(BW2Blue, 'holes');
BWnobordBlue = imclearborder(BWdfillBlue, 4);
seD = strel('diamond',1);
BWfinalBlue = imerode(BWnobordBlue,seD);
BWfinalBlue = imerode(BWfinalBlue,seD);
BWoutline = bwperim(BWfinalBlue);
SegoutBlue = rgbBlue;
SegoutBlue(BWoutline) = 255;
 
figure(1);
subplot(2,2,1),imshow(BWfinalBlue);
title('Subplot 1: Inner Area (Nucleus)')
subplot(2,2,2),imshow(BW1);
s1 = num2str(ratio);
s2 ='Index Nucleus/Citoplasm Analysis Ratio %';
s = strcat(s2,s1);
text(500, 2000, s)
title('Subplot 2: Outer Area (Citoplasm)')
 
%%[xmin ymin width height]
ZoomNucleus = imcrop(BWfinalBlue,[531 1681 320 240]);
ZoomCitoplasm = imcrop(BW1,[531 1681 320 240]);
 
subplot(2,2,3),imshow(ZoomNucleus);
title('ZOOM AREA (Nucleus)')
subplot(2,2,4),imshow(ZoomCitoplasm);
title('ZOOM AREA (Citoplasm)')
%%% mostrando 4 cuadrantes con imagne subrayada
figure(2);
subplot(2,2,1),imshow(rgb);
title('1. Original Image')
subplot(2,2,2),imshow(Segout);
title('2. Citoplasm and Nucleus delimited')
 
%%[xmin ymin width height]
ZoomDelimited = imcrop(Segout,[531 1681 320 240]);
ZoomGreenout = imcrop(rgbgreenout,[531 1681 320 240]);
 
subplot(2,2,3),imshow(ZoomDelimited);
title('3. ZOOM DELIMITED AREA')
subplot(2,2,4),imshow(ZoomGreenout);
title('4. ZOOM DELIMITED AREA HIGH CONTRAST')
%%ORINTATION CELLS DETECT
figure, imshow(SegoutBlue), title('Cell Polarity (Orientation Angle)');
 
[LBlue,NeBlue]=bwlabel(BWfinalBlue);
propied = regionprops(LBlue,'all');
hold on
%%for n=1:size(propied,1)
%%rectangle('position',propied(n).BoundingBox,'EdgeColor','g','LineWidth',2);
%%end
t = linspace(0,2*pi,50);
s=find([propied.Area]>=500&[propied.Area]<=5000);
for n=1:size(s,2)
%%rectangle('position',propied(s(n)).BoundingBox,'EdgeColor','r','LineWidth',2);
    a = propied(s(n)).MajorAxisLength/2;
    b = propied(s(n)).MinorAxisLength/2;
    Xc = propied(s(n)).Centroid(1);
    Yc = propied(s(n)).Centroid(2);
    phi = (-propied(s(n)).Orientation * pi)/180;%%degree to radians
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'r','Linewidth',1)
end
