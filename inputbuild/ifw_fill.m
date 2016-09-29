function [FRESHWATER] = ifw_fill(melt_total,month,new_mask,layer_vec)

   
    total_cells = sum(sum(new_mask));
    
    month_vec = [31 28 31 30 31 30 31 31 30 31 30 31];
    days = month_vec(month);
    
    fw_per_cell = melt_total/total_cells;
    
    FRESHWATER = zeros(60,384,320);
    indices = find(new_mask==1);
    [row,col] = ind2sub(size(new_mask),indices);
    z_val = find(layer_vec(:,2)~=0);
    
    for nn=1:length(indices);
        current_z = z_val(nn);
        FRESHWATER(layer_vec(current_z,2),row(nn),col(nn))=1;
        
    end
    
    FRESHWATER = FRESHWATER.*fw_per_cell;
    

end