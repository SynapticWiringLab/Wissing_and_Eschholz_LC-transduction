function [AvBright1, AvBright2, co_express] = f03_coexpress(CP_data, mastermask, overlap, im_template, im_second)


co_express = zeros(1, size(CP_data, 1));

for cellIdx = 1:size(CP_data, 1)
    cur_data = CP_data(cellIdx, :);
    x_data = cur_data(1:2:end); x_data = x_data(~isnan(x_data)); x_data = [x_data x_data(1)];
    y_data = cur_data(2:2:end); y_data = y_data(~isnan(y_data)); y_data = [y_data y_data(1)];
    cur_mask = poly2mask(x_data,y_data,size(im_template, 1),size(im_template, 2));
    
    AvBright1(cellIdx) =  mean(im_template(cur_mask));
    AvBright2(cellIdx) =  mean(im_second(cur_mask));
    
    if sum(sum(mastermask(cur_mask))) >= sum(sum(cur_mask))*overlap
        co_express(cellIdx) = 1;
    end
    
end


end

