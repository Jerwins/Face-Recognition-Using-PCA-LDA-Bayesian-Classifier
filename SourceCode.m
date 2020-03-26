
function [aimg alex aley arex arey amx amy] =imalign (bigimg, imheight, imwidth, filesavepath, leyex, reyex, eyey, mouthy)
if (nargin<3)
    Imheight=292;
    Imwidth=240;
    Leyex=56;
    Reyex=184;
    Eyey=88;
    Mouthy=224;
End
If nargin>2 && nargin<7
    Leyex=imwidth/4;
    Reyex=imwidth*3/4;
    Eyey=imheight/3;
    Mouthy=imheight*3/4;
End
Img=imread (bigimg);
Mouthx= (leyex+reyex)/2;
Os=size(img);
Figure; imshow (img);
[x,y]=getpts;
% align image
%compute rot scale deltX deltY
    lex=x (1);
    ley=y (1);
    rex=x(2);
    rey=y(2);
    mx=x(3);
    my=y(3); 
    hold on;
    plot([lex rex mx],[ley rey my],'bh');  
    X=[lex rex mx];
    Y=[ley rey my];
    Xc=[leyex reyex mouthx];
    Yc=[eyey eyey mouthy];
    A1=X*X'+Y*Y';
    A2=lex+rex+mx;
A3=ley+rey+my;
A4=X*Xc'+Y*Yc';
    A5=X*Yc'-Y*Xc';
    A6=leyex+reyex+mouthx;
    A7=eyey+eyey+mouthy;
    %At=b, t=A\b
    A=[A1 0 A2 A3;0 A1 -A3 A2;A2 -A3 3 0;A3 A2 0 3];
    b=[A4 A5 A6 A7]';
    t=A\b;
    k1=t(1,1);
    k2=t(2,1);
    t1=t(3,1);
    t2=t(4,1);
   rot=atan2(k2,k1);
    scale=(k1*k1+k2*k2)^0.5;
    ox=-t1;
    oy=-t2;
   angle=-rot*180/pi;
    aimg=imrotate(img,angle,'bilinear');
    % after rotation, the (0,0) is changed; 
    % calculate the shifts to adjust the crop
    sx=0;sy=0;
    if 0<angle && angle<=90 sy=os(2)*sin(-rot); end
    if -90<=angle && angle<0 sx=os(1)*sin(rot); end
    if 90<=angle && angle<=180 
        sx=os(2)*sin(-rot-pi/2); 
        sy=size(aimg,1);
    end
    if -90>=angle && angle>=-180 
        sx=size(aimg,2);
        sy=os(1)*sin(rot-pi/2);
    end  
    sx=sx*scale;
    sy=sy*scale;
    aimg=imresize(aimg,scale,'bilinear');
        cx=round(ox+sx+1);
    cy=round(oy+sy+1);
    aimg=imcrop(aimg,[cx,cy,imwidth-1,imheight-1]); 
    %if img not big enough, pad zeros
    if size(aimg,1)<imheight || size(aimg,2)<imwidth
        timg =cast (zeros (imheight, imwidth, size (aimg, 3)),class(aimg));
        if cx<1 px=-cx+1; else px=1; end
        if cy<1 py=-cy+1; else py=1; end
        timg(py:py+size(aimg,1)-1,px:px+size(aimg,2)-1,:)=aimg;
        aimg=timg;
end    
    alex=(lex*cos(rot)-ley*sin(rot))*scale-ox;
    aley=(lex*sin(rot)+ley*cos(rot))*scale-oy;
    arex=(rex*cos(rot)-rey*sin(rot))*scale-ox;
    arey=(rex*sin(rot)+rey*cos(rot))*scale-oy;
    amx=(mx*cos(rot)-my*sin(rot))*scale-ox;
    amy=(mx*sin(rot)+my*cos(rot))*scale-oy;
   figure;
   imshow(aimg)
    aimg=rgb2gray(aimg);
   colormap(gray(256));
 imwrite(aimg,filesavepath);
 hold on;
%plot ([alex arex amx],[aley arey amy],'r+');
    IntensityDif(filesavepath);
function ans= IntensityDif (inputimage)

M=60; 
filepath='Training set\';
indexa='a'; 
indexb='b'; 
img_matrix=[]; 
ImgDifVecs=[]; 
Vectors=[]; 
Values=[]; % To hold the eigenvalues (sorted) from pc_evectors
Psi=[];% To hold the mean image difference vector from pc_evectors
IntraCovMat=[]; % To hold the covariance matrix of the intrapersonal difference vectors
DimDesired =25; % The number of top eigenvectors I want
DiagEigValMat=[]; % To store the diagonal matrix of the top DimDesired Eigenvalues
imgMatrixA=[]; %To store 'a' filename images
imgMatrixB=[];%To store 'b' filename images
ExtraDifImgMatrix=[];%To hold the extrapersonal image difference vectors
ExtraDesiredSubDim=25; % The number of top extrapersonal eigenvalues 
for i=1:M
    strimgnameA=strcat(int2str(i), indexa);
    strfilenamea=strcat(strimgnameA, '.png');
    strfilepath=strcat(filepath, strfilenamea);
    img=imread(strfilepath);
    sizeA=numel(img);
    img=reshape(img', sizeA,1);
    img_matrix =[img_matrix img];
    imgMatrixA =[imgMatrixA img];
    strimgnameB=strcat(int2str(i), indexb);
    strfilenamea=strcat(strimgnameB, '.png');
    strfilepath=strcat(filepath, strfilenamea);
    img=imread(strfilepath);
    sizeA=numel(img);
    img=reshape(img', sizeA,1);
    img_matrix =[img_matrix img];
    imgMatrixB =[imgMatrixB img];
     end
for i =1:M
    ImgDifVecs(:,i) = imgMatrixA(:,i) - imgMatrixB(:,i);
ImgDifVecs(:,i+1) = imgMatrixB(:,i) - imgMatrixA(:,i);
   end
 [Vectors,Values,Psi,IntraCovMat] = pc_evectors(ImgDifVecs,DimDesired); 
i=30;
 [normalizing_denom] = cal_normalizingDenom(Values,DimDesired);
 [transVecMat]= perform_whitening(Vectors, Values,img_matrix,DimDesired);
imgtest= imread(inputimage);
sizeA=numel(imgtest);
imgtest=reshape(imgtest', sizeA,1);
[whitenedImgVec]= perform_whitening(Vectors, Values,imgtest,DimDesired);
 [MLresult,index]= cal_Likelihoods (transVecMat,whitenedImgVec,normalizing_denom);
if MLresult>1.3119e-063
    fprintf('recognized')
    msgbox('RECOGNIZED', 'Output')
else
    fprintf('not recognized');
    msgbox('Not recognized', 'Output')
    return
end
if(mod(index,2)==0)
    strimgname=strcat(int2str(index/2),'b');
    else
    strimgname=strcat(int2str(index/2),'a');
    end
strimgname
%Show the input image for feedback purpose
figure(3)
    strfilepath=inputimage;
img=imread(strfilepath);
    img=histeq(img,255);
subplot(ceil(sqrt(M)),ceil(sqrt(M)),i)
title('grayscale image')
    imshow(img)
drawnow;
%Show the image found for feedback purpose
figure(4);
strfilenamea=strcat(strimgname, '.png');
    strfilepath=strcat(filepath, strfilenamea);
    img=imread(strfilepath);
    img=histeq(img,255);
subplot(ceil(sqrt(M)),ceil(sqrt(M)),i)
title('recognized image')
    imshow(img)
   drawnow;
function [transVecMat]= perform_whitening(Vectors, Values, imgM,subDim)
DiagVal=[];
TempVal=[];
w=[];
for i =1:subDim
TempVal(i) = (Values(i))^(-1/2);
end
DiagVal= diag(TempVal);
transVecMat=[];%The matrix to hold the whitened image vectors
for i=1:size(imgM,2)
    transVecMat(:,i) = DiagVal * Vectors'* double(imgM(:,i)); %Perform whitening
end
function [normalizing_denom] = cal_normalizingDenom(Values,DimDesired) %L is the covariance matrix returned from pc_evectors
normalizing_const = 2*pi;
normalizing_const= normalizing_const ^(DimDesired/2); % This is the constant as specified in the first part of the normalizing formula
subVal=[];
 subVal=Values(1:DimDesired);
temp=1;
 for i =1:DimDesired
     temp =temp * subVal(i);
 end
normalizing_denom = normalizing_const * (temp)^(1/2)
function [normalizedWhitVecs] = cal_NormalizeWhit(normalizing_denom,transVecMat)
normalizedWhitVecs=[];
for i=1:size(transVecMat,2)
normalizedWhitVecs(:,i) = transVecMat(:,i)/normalizing_denom;
    end
function [MLresult,index]= cal_Likelihoods (WhitenVecs,ImgColVec,normalizing_denom)
likelihoods=[]
DBsize=size(WhitenVecs,2);
for i= 1:DBsize

   temp = ImgColVec - WhitenVecs(:,i);
    tempnorm = (-1/2)*((norm(temp))^2);
    likelihoods(i) = exp(tempnorm)/normalizing_denom;
   end
[MLresult,index]=max(likelihoods)
function [Vectors,Values,Psi,L] = pc_evectors(A,numvecs)
if nargin ~= 2
    error('usage: pc_evectors(A,numvecs)');
  end;
 nexamp = size(A,2);
fprintf(1,'Computing average vector and vector differences from avg...\n');
  Psi = mean(A')';
for i = 1:nexamp
    A(:,i) =A(:,i) - Psi;
  end;
       ele=size(A,2);
fprintf(1,'Calculating L=A''A\n');
  L=A'*A;
fprintf(1,'Calculating eigenvectors of L...\n');
  [Vectors,Values] = eig(L);
  fprintf(1,'Sorting evectors/values...\n');
  [Vectors,Values] = sortem2(Vectors,Values);
fprintf(1,'Computing eigenvectors of the real covariance matrix..\n');
  Vectors = A*Vectors;
  Values = diag(Values);
  Values = Values / (nexamp-1);
num_good = 0;
  for i = 1:nexamp
    Vectors(:,i) = Vectors(:,i)/norm(Vectors(:,i));
    if Values(i) < 0.0000001
      Values(i) = 0;
      Vectors(:,i) = zeros(size(Vectors,1),1);
    else
      num_good = num_good + 1;
end;
  end;
 if (numvecs > num_good)
    fprintf(1,'Warning: numvecs is %d; only %d exist.\n',numvecs,num_good);
    numvecs = num_good;
  end;
  Vectors = Vectors(:,1:numvecs);
function [vectors values] = sortem2(vectors, values)
if nargin ~= 2
 error('Must specify vector matrix and diag value matrix')
end;
vals = max(values); %create a row vector containing only the eigenvalues
[svals inds] = sort(vals,'descend'); %sort the row vector and get the indicies
vectors = vectors(:,inds); %sort the vectors according to the indicies from sort
values = max(values(:,inds)); %sort the eigenvalues according to the indicies from sort
values = diag(values); %place the values into a diagonal matrix








