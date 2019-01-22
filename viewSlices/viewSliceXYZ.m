function [  ] = viewSliceXYZ( x, n1, n2, n3, fid, filename, save, location )

if exist('fid', 'var') ~= 1
    fid = 1;
end
if exist('filename', 'var') ~= 1
    filename = 'x';
end
if exist('save', 'var') ~= 1
    save = 0;
end
if exist('location', 'var') ~= 1
    location = './figures/';
end 

img =  squeeze(x(:, n2, :));
[m, n] = size(img);

if m > n
   img = img';
end

figure(fid)
clf
colormap(gray)
imagesc(img)
title({filename, sprintf('z = %d', n2)}, 'interpreter', 'none')
axis image
colorbar
caxis([0 255])

if m < n
    xlabel('x');
    ylabel('y');
else
    xlabel('y');
    ylabel('x');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square

if save
    imwrite(uint8(img), sprintf('%s/%s.tiff',location, filename),'WriteMode', 'append');
end

%%
img = squeeze(x(:, :, n3));
[m, n] = size(img);

if m > n
   img = img';
end

figure(fid+1)
clf
colormap(gray)
imagesc(img)
axis image
colorbar
caxis([0 255])
title({filename, sprintf('y = %d', n3)}, 'interpreter', 'none')

if m < n
    xlabel('x');
    ylabel('z');
else
    xlabel('z');
    ylabel('x');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square

if save
    imwrite(uint8(img), sprintf('%s/%s.tiff',location, filename),'WriteMode', 'append');
end

%%
img = squeeze(x(n1, :, :));
[m, n] = size(img);

if m > n
   img = img';
end

figure(fid+2)
clf
colormap(gray)
imagesc(img)
axis image
colorbar
caxis([0 255])
title(sprintf('x = %d', n1))
title({filename, sprintf('x = %d', n1)}, 'interpreter', 'none')
if m < n
    xlabel('y');
    ylabel('z');
else
    xlabel('z');
    ylabel('y');
end

% axis off
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square

if save
    imwrite(uint8(img), sprintf('%s/%s.tiff',location, filename),'WriteMode', 'append');
end

end

