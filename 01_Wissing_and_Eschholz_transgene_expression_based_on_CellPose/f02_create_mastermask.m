function [master_mask] = f02_create_mastermask(CP_data, masksize)

master_mask = zeros(masksize);
for cellIdx = 1:size(CP_data, 1)
    cur_data = CP_data(cellIdx, :);
    x_data = cur_data(1:2:end); x_data = x_data(~isnan(x_data)); x_data = [x_data x_data(1)];
    y_data = cur_data(2:2:end); y_data = y_data(~isnan(y_data)); y_data = [y_data y_data(1)];
    cur_mask = poly2mask(x_data,y_data,masksize(1),masksize(2));
    master_mask = master_mask + cur_mask;
end


end

