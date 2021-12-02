function CofA=cofactor(A)

% Check Input Argument
if isempty(A)
    error(message('cof:EmptyMatrix','The matrix A does not exist'));
end

% Algorithm
[r, c]=size(A);     % number of rows and columns of A
temp_A=A;           % set a temporary matrix to calculate the determinants

% Initialize cofactor element matrix if A is not symbolic
sympar=1;
if ~isa(A,'sym')
    M=zeros(r,c);
    sympar=0;    
end    

switch sympar
%     if matrix A is symbolic than use the in-house function to calculate
%     the determinant of the matrix
    case 1
        for i=1:r,
            for j=1:c,
                temp_A(i,:)=[];             % remove i-th row
                temp_A(:,j)=[];             % remove j-th column
                detA=determin(temp_A);      % determinant of temporary matrix using the in-house function
                M(i,j)=((-1)^(i+j))*detA;   % compute cofactor element
                temp_A=A;                   % reset elements of temporary matrix to input elemens
            end
        end
%     if matrix A is numeric than use the Matlab built-in function to 
%     calculate the determinant of the matrix          
    case 0
        for i=1:r,
            for j=1:c,
                temp_A(i,:)=[];             % remove i-th row
                temp_A(:,j)=[];             % remove j-th column
                detA=det(temp_A);           % determinant of temporary matrix using the matlab built-in function
                M(i,j)=((-1)^(i+j))*detA;   % compute cofactor element
                temp_A=A;                   % reset elements of temporary matrix to input elemens
            end
        end
end

CofA=M;

end