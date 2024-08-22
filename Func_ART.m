function [p_est,f_temp,f_final]=Func_ART(p,w,lambda,f_final,f_temp,image_size)
for i=1:length(p)
    
        if sum(w(i,:))==0
            f_final=f_temp;
           
        else   
            qi(i)=(w(i,:)*f_temp);
            f_final=f_temp+lambda*((p(i)-qi(i))/sum(w(i,:).^2))*w(i,:)';
        
        end
        
        f_temp=f_final;
        p_est=reshape(f_final',[image_size,image_size])';
     end
%         disp([num2str(ii),' iteration done']);
end
        