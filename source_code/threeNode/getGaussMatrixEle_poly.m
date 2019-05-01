function [result] = getGaussMatrixEle_poly(i,j,x_0,h, p, q, scale)

c = 0.8;
phi = @(x)(exp(-c*x*x));
phi_d = @(x)(exp(-c*x*x)*(-2*c*x));



if mod(i,2)==0 && mod(j,2)==0  %% both are node

    
    if abs(i-j)> 2
        result = 0;
    elseif i == j && i == scale*2
        f1 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(i/2-1)*h,h), p);
        f2 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(i/2-1)*h,h), f1);
        
        result = integral(f2,(i/2-1)*h,i/2*h,'ArrayValued',true); 
        
        f3 = mulFunc(getShape_poly(phi,phi,3,(i/2-1)*h,h), q); 
        f4 = mulFunc(getShape_poly(phi,phi,3,(i/2-1)*h,h), f3); 
        result = result + integral(f4,(i/2-1)*h,i/2*h,'ArrayValued',true);
        
    elseif abs(i-j) == 2
        index = min(i,j);
        
        f1 = mulFunc(getShapeDeri_poly(phi,phi_d,1,(index/2)*h,h), p);
        f2 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(index/2)*h,h), f1);
        
        result = integral(f2,index/2*h,(index/2+1)*h,'ArrayValued',true); 
    
        f3 = mulFunc(getShape_poly(phi,phi,1,(index/2)*h,h), q); 
        f4 = mulFunc(getShape_poly(phi,phi,3,(index/2)*h,h), f3); 
        
        result = result + integral(f4,index/2*h,(index/2+1)*h,'ArrayValued',true);
        
    else
        f1 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(i/2)*h-h,h), p);
        f2 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(i/2)*h-h,h), f1);
        
        result = integral(f2,i/2*h-h,(i/2)*h,'ArrayValued',true); 
    
        f3 = mulFunc(getShape_poly(phi,phi,3,(i/2)*h-h,h), q); 
        f4 = mulFunc(getShape_poly(phi,phi,3,(i/2)*h-h,h), f3); 
        
        result = result + integral(f4,i/2*h-h,(i/2)*h,'ArrayValued',true);
        
        f5 = mulFunc(getShapeDeri_poly(phi,phi_d,1,(i/2)*h,h), p);
        f6 = mulFunc(getShapeDeri_poly(phi,phi_d,1,(i/2)*h,h), f5);
        
        result = result + integral(f6,i/2*h,(i/2+1)*h,'ArrayValued',true); 
    
        f7 = mulFunc(getShape_poly(phi,phi,1,(i/2)*h,h), q); 
        f8 = mulFunc(getShape_poly(phi,phi,1,(i/2)*h,h), f7); 
        
        result = result + integral(f8,i/2*h,(i/2+1)*h,'ArrayValued',true);
        
    end


elseif mod(i,2)==1 && mod(j,2)==1 %% both are mid
    if i == j
        f1 = mulFunc(getShapeDeri_poly(phi,phi_d,2,(i-1)*(h/2),h), p);
        f2 = mulFunc(getShapeDeri_poly(phi,phi_d,2,(i-1)*(h/2),h), f1);
        result = integral(f2, (i-1)*(h/2), (i-1)*(h/2)+h, 'ArrayValued',true);
        f3 = mulFunc(getShape_poly(phi,phi,2,(i-1)*(h/2),h), q);
        f4 = mulFunc(getShape_poly(phi,phi,2,(i-1)*(h/2),h), f3);
        result = result + integral(f4, (i-1)*(h/2), (i-1)*(h/2)+h, 'ArrayValued',true);
    else
        result = 0;
    end
else                              %% one mid, one node    
    if  j - i == 1
        if mod(i,2) == 0         % i, j = 2 3
            f1 = mulFunc(getShapeDeri_poly(phi,phi_d,1,(i)*(h/2),h), p);
            f2 = mulFunc(getShapeDeri_poly(phi,phi_d,2,(i)*(h/2),h), f1);
            result = integral(f2, (i)*(h/2), (i)*(h/2)+h, 'ArrayValued',true);
            f3 = mulFunc(getShape_poly(phi,phi,1,(i)*(h/2),h), q);
            f4 = mulFunc(getShape_poly(phi,phi,2,(i)*(h/2),h), f3);
            result = result + integral(f4, (i)*(h/2), (i)*(h/2)+h, 'ArrayValued',true);
        else                     % i, j = 1 2
            f1 = mulFunc(getShapeDeri_poly(phi,phi_d,3,(i)*(h/2)-h/2,h), p);
            f2 = mulFunc(getShapeDeri_poly(phi,phi_d,2,(i)*(h/2)-h/2,h), f1);
            result = integral(f2, (i)*(h/2)-h/2, (i)*(h/2)+h/2, 'ArrayValued',true);
            f3 = mulFunc(getShape_poly(phi,phi,3,(i)*(h/2)-h/2,h), q);
            f4 = mulFunc(getShape_poly(phi,phi,2,(i)*(h/2)-h/2,h), f3);
            result = result + integral(f4, (i)*(h/2)-h/2, (i)*(h/2)+h/2, 'ArrayValued',true);
        end
        
    elseif i - j == 1
        result = getGaussMatrixEle_poly(j,i,x_0,h, p, q, scale);                                     % remember to modifiy
    else
        result = 0;
    end
end