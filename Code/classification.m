%% Everything I wrote above in http://uk.mathworks.com/matlabcentral/answers/263708-hyperspectral-image-classification-unmixing-matlab-code#answer_207146 applies to hyperspectral data.

%% I did mention JPG files, but keep in mind that one possible...
%  representation of hyperspectral data is to have one image file per data layer.

%% Look at this code:
%%
%train the classifier
classification_information = train_classifier('SomeFileOfGroundTruth.xls');
%read the hyperspectral data
YourData = Hyperspectral_read('YourImageFile.dat');
[rows, cols, bands] = size(YourData);
%initialize the output matrix
class_matrix = zeros(rows, cols, 'uint8');  %use 'uint16' if you have more than 255 classes of vegetation
%classify every pixel
for R = 1 : rows
  for C = 1 : cols
    this_pixel = squeeze(YourData(R, C, :));
    this_class_number = classify_pixel(this_pixel, classification_information);
    class_matrix(R, C) = this_class_number;
  end
end
%now create a colormap. In real code this would not be
%done dynamically: you would want to have particular colors
%for particular vegetation.
maxclass = max(class_matrix(:));
cmap = jet(maxclass);
cmap(1,:) = 0;         %black for no vegetation
%now display the class matrix
image(class_matrix);
colormap(cmap)
colorbar()
%% Is there anything in that code that requires that the data...
%  be an RGB matrix? At no time does it involve an RGB matrix.
%  The display portion does involve an indexed (pseudocolor) image,
%  each element of which is the class number of the vegetation as ...
%  determined by examining all of the bands with whatever information...
%  is returned by the training of the classifier.