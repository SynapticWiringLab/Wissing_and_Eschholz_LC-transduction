function  f01_drawoutline(cellpose_data, color)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for cellIdx = 1:size(cellpose_data, 1)
    plot(cellpose_data(cellIdx, 1:2:end), cellpose_data(cellIdx, 2:2:end), 'color', color./255, 'LineWidth', 1.5)
end

end

