%% ===================================
% Cálculo das celeridades

function mat = material(E,nu,rho,d)

mat.E = E;
mat.nu = nu;
mat.rho = rho;
mat.d = d;

mat.lambda 	= E*nu/((1+nu)*(1-2*nu));
mat.mu		= E/(2*(1+nu));

mat.cL 		= sqrt((mat.lambda+2*mat.mu)/rho);
mat.cT		= sqrt(mat.mu/rho);

end
%% ===================================