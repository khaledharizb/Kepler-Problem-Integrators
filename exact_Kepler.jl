function phiKepler3D(z,h) 
tol=1.e-14;  q₀ = z[1:3]; p₀ = z[4:6]; r₀ = sqrt(q₀' * q₀);  
H₀ = 0.5 * p₀' * p₀ - 1 / r₀ ; L₀ = cross(q₀,p₀) ; RL = cross(p₀,L₀) - q₀/r₀; 
a = -1 / (2*H₀) ; w = a^(3/2); ECC = sqrt(1 + 2 * H₀ * L₀' * L₀);
 
g(E) = ECC * sin(E) + h / w ;
  x₀ = w * h;  E = FixIter(g,x₀,tol);
C1 = a * RL / ECC ; C2 =  cross(sqrt(a) * L₀,RL/ ECC) ; 
    
q =  C1 * cos(E) +  C2 *  sin(E) - a * RL;
 dsdt = w * (1 - ECC * cos(E));    
p =  - C1 * sin(E) / dsdt  +  C2 *  cos(E) / dsdt;
       
 return [q;p];

end

function phiKepler3D(z,h,steps) 
q = zeros((size(z)[1],steps+1));
q[:,1] = z;    
t0  = 0;
for i=1:steps;
  t0 = t0+h;
 q[:,i+1] = phiKepler3D(z,t0) ; 
end   
    return q;
end;
