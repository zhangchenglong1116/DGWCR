% compute squared Euclidean distance
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
function d = Get_Distance(a)
    sm=ones(1,size(a,1));
    aa=sm*(a.*a);  ab=(a'*a); 
    aa = round(10000*aa)/10000;    ab = round(10000*ab)/10000;     
    d = repmat(aa',[1 size(aa,2)]) + repmat(aa,[size(aa,2) 1]) - 2*ab;
    d = real(d);d = sqrt(max(d,0));
end