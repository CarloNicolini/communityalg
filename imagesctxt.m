function imagesctxt(mat,clim)
if nargin==1
    clim=[];
end

%IMAGESCTXT Scale data and display as image with text of values impressed.
%   IMAGESCTXT(...) is the same as IMAGESC(...) except the values are
%   printed on pixels.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

imagesc(mat,clim);
colormap(parula);

textStrings = num2str(mat(:),'%1.3g');  % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x,y] = meshgrid(1:size(mat,2),1:size(mat,1));   % Create x and y coordinates for the strings

% Plot the strings
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  % Get the middle value of the color range

% text color of the strings so they can be easily seen over the background color
textColors = repmat(mat(:) < midValue,1,3);

set(hStrings,{'Color'},num2cell(textColors,2));  % Change the text colors

set(gca,'XTick',1:size(mat,2),...                         % Change the axes tick marks
        'XTickLabel',num2cell(1:size(mat,2)),...  %   and tick labels
        'YTick',1:size(mat,1),...
        'YTickLabel',num2cell(1:size(mat,1)),...
        'TickLength',[0 0]);