function [Fx, Fy] = Get_Comp_Fields(omega, BC, Dim,h, Fz)
    N = round(Dim(1)/h);%num of x dim grid points
    M = round(Dim(2)/h);%num of y dim grid points
    
    temp = Fz.';
    Fz = temp(:);
    
    switch BC{1}{1}
        case 'pml'
            A = -ones(M*N,2)/h;
            A(:,2) = 1/h;
            %skip first element
            %positive offset skip first few elements
            %negative offset skip last few elements
            A(N+1:N:(M*N+1),2) = 0; 
            Dxe = spdiags(A,[0 1], M*N,M*N);
        case 'periodic'
            A = -ones(M*N,3)/h;
            A(:,3) = 1/h;
            A(N+1:N:(M*N+1),3) = 0;
            A(:,1) = 0;
            A(1:N:(M*N+1),1) = 1/h;
            Dxe = spdiags(A,[-(N-1) 0 1], M*N,M*N);
        otherwise
            warning('invalid x boundary condition');
    end

    switch BC{2}{1}
        case 'pml'
            A = -ones(M*N,2)/h;
            A(:,2) = 1/h;
            Dye = spdiags(A,[0 N], M*N,M*N);
        case 'periodic'
            A = ones(M*N,3)/h;
            A(:,2) = -1/h;
            Dye = spdiags(A,[N-M*N 0 N], M*N,M*N);
        otherwise
            warning('invalid y boundary condition');
    end
    
    %Ce = [Dye; -Dxe];
    Fx = (1./(1i*omega)).*(Dye*Fz);
    Fy = -(1./(1i*omega)).*(Dxe*Fz);
    
    Fx = reshape(Fx,[N,M]).';
    Fy = reshape(Fy,[N,M]).';
end