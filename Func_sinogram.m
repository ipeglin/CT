function [projection_data_new]=Func_sinogram(projection_data,image_size)
if image_size==128;
    ray_number=size(projection_data,1);
    projection_data_new=projection_data((ray_number-183)/2+1:(ray_number-(ray_number-183)/2),:);
else if image_size==50
        ray_number=size(projection_data,1);
        projection_data_new=projection_data((ray_number-73)/2+1:(ray_number-(ray_number-73)/2),:);
    end

end
end