function I_CoF = CoOccurrenceFilter()

%Implementation of Co-Occurrence Filter
I = imread('TestImage.png');
I_CoF = zeros(size(I), 'uint8');

%% Parameters
filter_size = 15;
sigma = 2*(filter_size)^0.5 + 1; %For debug: 1e-12 , should return original image

%Calculate normalized co-occurrence matrix
M = collectOccurrenceStatistic(I,filter_size,sigma);

%Create gaussian mask
h = createGaussianFilter(filter_size,sigma);

for color= 0:255
    
    [idx_y, idx_x] = find(I==color);
    
    for k=1:length(idx_x)
        
        %Extract window center
        window_center = [idx_x(k) idx_y(k)]; %In format: [x y]
        
        %Cut patch from image and the the gaussian mask to its size
        [patch, h_current] = fitMaskToIdx(double(I),h,window_center,filter_size);
        
        %eq. (5) from the paper
        w = h_current.*reshape( M(color+1, patch(:)+1 ) , size(patch));
        
        %eq. (1) from the paper
        I_CoF(window_center(2),window_center(1)) =  uint8(sum2(w.*double(patch))/sum2(w));
        
    end
end

%For debug
figure; 
subplot(1,2,1); imshow(I); title('I') ; 
subplot(1,2,2); imshow(I_CoF); title('I CoF');

end



function [M] = collectOccurrenceStatistic(I,filter_size,sigma)
%Inputs: I- gray image on which the statistic will be collected
%        filter_size - odd number, window size for which the statistic will be collected for each pixel
%        sigma- variance of the distance gaussian
%Output: M- normalized co-occurrence matrix. eq. (6) from the article

%Normalization factor
I_histogram = imhist(uint8(I));
h_a_h_b = I_histogram*I_histogram';


%Create gaussian mask. histogram weight
h = createGaussianFilter(filter_size,sigma);

%Initialize co-occurrence matrix
C = zeros(256,256);
for color = 0:255
    [idx_y, idx_x] = find(I==color);
    
    %Calculate the co-accurance in a specific window
    c_per_level = zeros(256,1);
    for k=1:length(idx_x)
        
        %Extract window center
        window_center = [idx_x(k) idx_y(k)]; %In format: [x y]
        
        %Cut patch from image and the the gaussian mask to its size
        [patch, h_current] = fitMaskToIdx(double(I),h,window_center,filter_size);
        
        %Accumulate- single elemnt in eq (7) in the paper
        c_per_level = c_per_level + accumarray(patch(:)+1,h_current(:), [256 1]);
        
    end
    
    C(color+1, :) = c_per_level;
end


%Normalization. eq (6) in the paper
M = C ./ (h_a_h_b + eps);

end


function [patch,h_current] = fitMaskToIdx(I, h, window_center, filter_size)
%Extract patch from the image based on the center and window size. fit also the mask
%Inputs:    I- image from which the patch will be extracted
%           h - mask which need to be fitted to the patch
%           window_center- center of the patch
%           filter_size- size of the patch.
%Outputs:   patch- cropped image with the size [filter_size filter_size ]
%           h_currnet - cropped mask which fits to the patch
image_size = size(I);
window_center_x = window_center(1);
window_center_y = window_center(2);

%calculate minimum and maxmimum in order to cut the image into patch
patch_idx_min_x = ceil(window_center_x - filter_size/2);
patch_idx_max_x = floor(window_center_x + filter_size/2);
patch_idx_min_y = ceil(window_center_y - filter_size/2);
patch_idx_max_y = floor(window_center_y + filter_size/2);


%If patch is outside the image- cut it along with the gaussian
mask_min_x = 0; mask_min_y = 0; mask_max_x = 0; mask_max_y = 0;
if patch_idx_min_x<1
    mask_min_x = 1 - patch_idx_min_x;
    patch_idx_min_x = 1;
end
if patch_idx_min_y<1
    mask_min_y = 1 -patch_idx_min_y;
    patch_idx_min_y = 1;
end
if patch_idx_max_x>image_size(2)
    mask_max_x = patch_idx_max_x - image_size(2);
    patch_idx_max_x = image_size(2);
end
if patch_idx_max_y>image_size(1)
    mask_max_y =  patch_idx_max_y - image_size(1);
    patch_idx_max_y = image_size(1);
end

%Assign to output
patch = I(patch_idx_min_y:patch_idx_max_y,patch_idx_min_x:patch_idx_max_x);
h_current = h(1+mask_min_y:end-mask_max_y,1+mask_min_x:end-mask_max_x);

end

function [h] = createGaussianFilter(filter_size,sigma)
%Creating gussian mask
%Inputs:    filter_size- mask size
%           sigma- variance of the gaussian
%Outputs:   h- array with vlues of the gaussian of the size [filter_size filter_size]
[X, Y] = meshgrid(ceil(-filter_size/2):floor(filter_size/2),ceil(-filter_size/2):floor(filter_size/2));
h = exp(-(X.^2 + Y.^2) / (2*sigma^2+eps));
end