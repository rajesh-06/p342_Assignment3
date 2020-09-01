#introducing partial pivoting to interchange 2 row if diagonal elements is zero
def par_pivot(A,B):
	n=len(A)
	for r in range(n):
		if A[r][r]==0:
			for r1 in range(r+1,n):
				if abs(A[r1][r])>A[r][r] and A[r][r]==0:
					(A[r],A[r1])=(A[r1],A[r])
					(B[r],B[r1])=(B[r1],B[r])
				else:
					continue
			else:
				continue
#Gauss-Jordan elimination to get solutions of linear equations
def gauss(A,B):
	m=len(A)
	n=len(A[0])
	for r in range(m):
		par_pivot(A,B)
		pivot=A[r][r]
		for c in range(r,n):
			A[r][c]=A[r][c]/pivot
		B[r]=B[r]/pivot
		for r1 in range(m):
			if r1==r or A[r1][r]==0:
				continue
			else:
				factor=A[r1][r]
				for c in range(r,n):
					A[r1][c]=A[r1][c]-A[r][c]*factor
				B[r1]=B[r1]-B[r]*factor
				
				
				
#matrix multiplication of 3*3 and 3*1
def Mat_multiplication(A,B):
	p=[0,0,0]
	for i in range(3):
		for j in range(3):
			p[i] = p[i] + (A[i][j] * B[j])	
	print("Ax =",p)
#matrix multiplication of 3*3 and 3*3
def Sqr_Mat_mult(A,B):
	p1=[[0,0,0], [0,0,0], [0,0,0]]	
	for i in range(3):
		for x in range(3):
			for j in range(3):
				p1[i][j] = p1[i][j] + (A[i][x]*B[x][j]) 
	return(p1)
#-------------------------------
#Solution for Q1
#importing matrix from txt file
print('Solution for Q1:')
with open('Q1_A.txt', 'r') as f:
    A = [[float(num) for num in line.split(',')] for line in f]
print('A = \n', A)

with open('Q1_B.txt', 'r') as f:
	for line in f:
		B=[float(num) for num in line.split(',')]
	print('B =', B)

#finding the solution
gauss(A,B)
print('x =', B)

#checking product
with open('Q1_A.txt', 'r') as f:
    A1 = [[float(num) for num in line.split(',')]for line in f]
Mat_multiplication(A1,B)
print('\n')


#-------------------------------
#solution for Q2
print('Solution for Q2:')
#importing matrix from txt file
with open('Q2_A.txt', 'r') as f:
	C =[[float(num) for num in line.split(',')]for line in f]
print('A = \n', C)

with open('Q2_B.txt', 'r') as f:
    for line in f:
    	D =[float(num) for num in line.split(',')]
print('B =', D)

#finding the solution
gauss(C,D)
print('x =', D)

#checking product
with open('Q2_A.txt', 'r') as f:
    C1 = [[float(num) for num in line.split(',')]for line in f]
Mat_multiplication(C1,D)
print('\n')

#-------------------------------
#solution for Q3
#importing matrix from txt file
print('Solution for Q3:')
with open('Q3_A.txt', 'r') as f:
	e =[[float(num) for num in line.split(',')]for line in f]
print('A = \n', e)
#Storing(deepcopy) the exact value of A
import copy
A1 = copy.deepcopy(e)
A2 = copy.deepcopy(e)
A3 = copy.deepcopy(e)
A4 = copy.deepcopy(e)
I=[[1,0,0],[0,1,0],[0,0,1]]
#calling Gauss-Jordan to transform the individual rows of identity matrix
gauss(A1,I[0])
gauss(A2,I[1])
gauss(A3,I[2])
#As transformed(I) is the transpose of A^(-1) we need to transpose it 
for i in range(3):
	for j in range(3):
		A4[i][j] = I[j][i]
#printing A^(-1)
print('A^(-1) =',A4)
#checking A*A^(-1) = I
print('A*A^(-1) =',Sqr_Mat_mult(e,A4))
#checking A^(-1)*A = I
print('A^(-1)*A =',Sqr_Mat_mult(A4,e))
#-------------------------------

'''Output of the code:


Solution for Q1:
A =
 [[1.0, 3.0, 2.0], [2.0, 7.0, 7.0], [2.0, 5.0, 2.0]]
B = [2.0, -1.0, 7.0]
x = [3.0, 1.0, -2.0]
Ax = [2.0, -1.0, 7.0]


Solution for Q2:
A =
 [[0.0, 2.0, 5.0], [3.0, -1.0, 2.0], [1.0, -1.0, 3.0]]
B = [1.0, -2.0, 3.0]
x = [-2.0, -2.0, 1.0]
Ax = [1.0, -2.0, 3.0]


Solution for Q3:
A =
 [[1.0, -3.0, 7.0], [-1.0, 4.0, -7.0], [-1.0, 3.0, -6.0]]
A^(-1) = [[-3.0, 3.0, -7.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]]
A*A^(-1) = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
A^(-1)*A = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

[Program finished]'''