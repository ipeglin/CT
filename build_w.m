function [weighting_factor]=build_w(w,projection_number,image_size,ray_number)

weighting_factor=zeros(projection_number*ray_number,image_size^2);
delta_theta=180/projection_number;
for k=1:projection_number
    weighting_factor(ray_number*(k-1)+1:ray_number*k,:)=w(ray_number*(delta_theta*k-1)+1:ray_number*(delta_theta*k),:);

end