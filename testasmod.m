clear all;

A=zeros(6,6);

A(1,2)=1;
A(1,3)=1;
A(1,5)=1;
A(2,3)=1;
A(3,4)=1;
A(4,5)=1;
A(5,6)=1;
A(2,4)=1;

A=triu(A);
A=A+A';
dlmwrite('test.adj',A,' ');

fprintf('--OLD--\n');
qold = modularity(A,[1 1 1 2 2 2])
fprintf('--NEW--\n');
qnew = modularity(A,[1 1 1 1 2 2])

