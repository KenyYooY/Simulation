//============================================================================
// Name        : Introduction.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define _DEBUG_ 3

#include <iostream>
#include <cmath>
#include <cstdlib>
//#include "Polynominal.h"

using namespace std;

void Print_Array(double *array,int n){

	for(int i=0; i<n; i++){
			cout << array[i] << ",";
	}
	cout << endl;

}

void Print_Array_2D(double *array,int row, int col){

	for(int i=0; i<row; i++){
		for(int j=0; j<col; j++){
			cout << array[i*col+j];
			if(j<col) cout << ",";
		}
		cout << endl;
	}

}

double Check_DiagonalPriority(double *a, int n){

	double Mi=0.0;

	for(int i=0; i<n; i++){
		double M=0.0;
		for(int j=0; j<n; j++){
			if(i!=j) M += abs(a[i*n+j]);
		}
		M /= abs(a[i*n+i]);
		if(M>Mi) Mi=M;
	}

	return(Mi);

}


int Pivot(double *a, double *b, int n){

	double adia_max;
	double err=1e-6;

	for(int i=0; i<n; i++){
		int i_max=i;
		adia_max=abs(a[i*n+i]);

		if(i+1<n){
			for(int ii=i+1;ii<n; ii++){
				double adia_max_tmp = abs(a[ii*n+i]);
				if(adia_max_tmp>adia_max){
					i_max=ii;
					adia_max=adia_max_tmp;
				}
			}
		}
		// Error : Too Small Diagonal
		if(adia_max<=err){
			cout << "K=" << i << ", FMAX=" << adia_max << endl;
			return -1;
		}
		// Exchange :
		if(i_max!=i){
			for(int j=i; j<n; j++){
	            double aa=a[i_max*n+j];
	            a[i_max*n+j]=a[i*n+j];
	            a[i*n+j]=aa;
			}
            double bb=b[i_max];
            b[i_max] = b[i];
            b[i] = bb;
		}
	}

	return 1;

}

/*
 * 初等解法による連立１次方程式の解：ガウス・ジョルダンの消去法
 */
int GaussJordan(double *a, int n, double *x){

	double fmax;
	double e=0.000001;
	int n1=n+1;

	for(int i=0; i<n; i++){
		// Pivot
		int ma=i;
		fmax=abs(a[i*n1+i]);

		if(i+1<n){
			int i1=i+1;
			for(int ii=i1;ii<n; ii++){
				double fmax1 = abs(a[ii*n1+i]);
				if(fmax1>fmax){
					ma=ii;
					fmax=fmax1;
				}
			}
		}
		// Error : Too Small Diagonal
		if(fmax<=e){
			cout << "K=" << i << ", FMAX=" << fmax << endl;
			return -1;
		}
		// Exchange :
		if(ma!=i){
			for(int j=i; j<n1; j++){
	            double aa=a[ma*n1+j];
	            a[ma*n1+j]=a[i*n1+j];
	            a[i*n1+j]=aa;
			}
		}
		// Solve
		double b=a[i*n1+i];

		for(int j=1; j<n1; j++){
            a[i*n1+j]=a[i*n1+j]/b;
		}
		for(int l=0; l<n; l++){
			if(l!=i){
				double s=a[l*n1+i];
				for(int j=0; j<n1; j++){
					a[l*n1+j]=a[l*n1+j]-s*a[i*n1+j];
				}
			}
		}

		for(int i=0; i<n; i++){
			x[i]=a[i*n1+(n1-1)];
		}
	}

	return 1;

}

void Test_GaussJordan(){
// Test by WolframAlpha : x+3y-5z=-9,3w-2x+y+z=6,8w-5x-7y+2z=-15,w+x-2y-3z=-15
	const int row = 4;
	double array[row][row+1]={
			{ 0.0 , 1.0 , 3.0 ,-5.0 ,-9.0},
			{ 3.0 ,-2.0 , 1.0 , 1.0 , 6.0},
			{ 8.0 ,-5.0 ,-7.0 , 2.0 ,-15.0},
			{ 1.0 , 1.0 ,-2.0 ,-3.0 ,-15.0}
	};

	double x[row];

	Print_Array_2D(array[0], row, row+1);

	if(GaussJordan(array[0], row, x)>0){
		Print_Array_2D(x,row,1);
	}else{
		cout << "Error" << endl;
	}

}

int InverseMatrixByGaussJordan(double *a,int n, double *ainv){
	// int n1=n+1;
	int m=n*2;
	double e=0.000001;
	double aainv[n][m];
	double fmax;

	// Copy component of matrix and set unit matrix into aainv
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			if(j<n){
				aainv[i][j]=a[i*n+j];
			}else if(i+n==j){
				aainv[i][j]=1.0;
			}else{
				aainv[i][j]=0.0;
			}
		}
	}
#if _DEBUG_
	Print_Array_2D(aainv[0], n, m);
#endif

	// Main Routin
	for(int i=0; i<n; i++){
		//// Pivot
		int ma=i;
		fmax=abs(aainv[i][i]);

		if(i+1<n){
			int i1=i+1;
			for(int ii=i1; ii<n;  ii++){
				double fmax1 = abs(aainv[ii][i]);
				if(fmax1>fmax){
					ma=ii;
					fmax=fmax1;
				}
			}
		}
		// Error : Too Small Diagonal
		if(fmax<=e){
			cout << "K=" << i << ", FMAX=" << fmax << endl;
			return -1;
		}
		// Exchange :
		if(ma!=i){
			for(int j=i; j<m; j++){
	            double aa=aainv[ma][j];
	            aainv[ma][j]=aainv[i][j];
	            aainv[i][j]=aa;
			}
		}

		// Solve
		double b=aainv[i][i];

		for(int j=i; j<m; j++){
			aainv[i][j]=aainv[i][j]/b;
		}
		for(int l=0; l<n; l++){
			if(l!=i){
				double s=aainv[l][i];
				for(int j=0; j<m; j++){
					aainv[l][j]=aainv[l][j]-s*aainv[i][j];
				}
			}
		}
	}


	// Copy inverse matrix to ainv from aainv
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			ainv[i*n+j]=aainv[i][j+n];
		}
	}
	return 1;
}

void Test_InverseMatrixByGaussJordan(){
// Test by WolframAlpha : inverse matrix[{0,2,3},{4,0,5},{1,2,0}]
	const int row = 3;
	double a[][row]={
			{ 0.0 , 2.0 , 3.0 },
			{ 4.0 , 0.0 , 5.0 },
			{ 1.0 , 2.0 , 0.0 }
	};

	double ainv[row][row];

	Print_Array_2D(a[0], row, row);

	if(InverseMatrixByGaussJordan(a[0], row,  ainv[0])>0){
		Print_Array_2D(ainv[0],row,row);
	}else{
		cout << "Error" << endl;
	}

}

/*
 * 初等解法による連立１次方程式の解：LU分解
 * 連立方程式Ax=b の係数行列Aを上三角行列U,下三角行列Uに分解(A=L・U)して
 * 連立方程式を L・Ux=b とし、 Ux=x'とおき、Lx'=b とする。
 * 下三角行列Lは既知のためLx'=bは一番上から解いていけて、x'が求まる。
 * 上三角行列Uも既知のためUx=x'は一番下から解いていけて、xが求まる。
 */
int LU_Decomposition(double *a, int n, double *b, double *x){

	double l[n][n], u[n][n];

	//// Calc. L,U Matrix
	// Column 1 (L[i,1], U[i,1])
	for(int i=0; i<n; i++){
		int j=0;
		// L[i,1]
		l[i][j] = a[i*n+j];
		// U[i,1]
		if(i>j){
			u[i][j] = 0.0;
		}else if(i==j){
			u[i][j] = 1.0;
		}
	}

	for(int j=1; j<n; j++){
		for(int i=0; i<n; i++){

			// Make U(Upper Triangle) Matrix
			if(i==0){ // make U[1,j]
				u[i][j] = a[i*n+j]/l[i][i];
			}else if(i==j){
				u[i][j] = 1.0;
			}else if(i>j){
				u[i][j] = 0.0;
			}else if(i<j){
				double sum_lu = 0.0;
				for(int k=0; k<i; k++) sum_lu += l[i][k]*u[k][j];
				u[i][j] = (a[i*n+j]-sum_lu)/l[i][i];
			}

			// Make L(Lower Triangle) Matrix
			if(i<j){
				l[i][j] = 0.0;
			}else{
				double sum_lu = 0.0;
				for(int k=0; k<j; k++) sum_lu += l[i][k]*u[k][j];
				l[i][j] = a[i*n+j]-sum_lu;
			}
		}
	}
#if _DEBUG_
	cout << "L Matrix is" << endl; Print_Array_2D(l[0], n, n);
	cout << "U Matrix is" << endl; Print_Array_2D(u[0], n, n);
#endif

	// Calc. x'

	double xp[n];
	xp[0]=b[0]/l[0][0];
	for(int i=1; i<n; i++){
		double sum_lxp = 0.0;
		for(int j=0; j<i; j++) sum_lxp += l[i][j]*xp[j];
		xp[i] = (b[i]-sum_lxp)/l[i][i];
	}
#if _DEBUG_
	cout << "Vector x' is" << endl; Print_Array_2D(xp, n, 1);
#endif

	// Calc. x
	x[n-1]=xp[n-1];
	for(int i=n-2; i>=0; i--){
		double sum_ux = 0.0;
		for(int j=i+1; j<n; j++) sum_ux += u[i][j]*x[j];
		x[i] = xp[i]-sum_ux;
	}

	return 1;
}

/*
 * L,U分解の行列をメモリ削減のため１つの行列にまとめた
 */
int LU_Decomposition2(double *a, int n, double *b, double *x){

	double lu[n][n];

	//// Calc. L,U Matrix
	// Column 1 (L[i,1], U[i,1])
	for(int i=0; i<n; i++){
		int j=0;
		// L[i,1]
		lu[i][j] = a[i*n+j];
	}

	for(int j=1; j<n; j++){
		for(int i=0; i<n; i++){

			// Make U(Upper Triangle) and L(Lower Triangle) Matrix
			if(i==0){ // make U[1,j]
				lu[i][j] = a[i*n+j]/lu[i][i];
			}else if(i<j){ // make U[i,j]
				double sum_lu = 0.0;
				for(int k=0; k<i; k++) sum_lu += lu[i][k]*lu[k][j];
				lu[i][j] = (a[i*n+j]-sum_lu)/lu[i][i];
			}else{ // make L[i,j]
				double sum_lu = 0.0;
				for(int k=0; k<j; k++) sum_lu += lu[i][k]*lu[k][j];
				lu[i][j] = a[i*n+j]-sum_lu;
			}
		}
	}
#if _DEBUG_
	cout << "L+U Matrix is" << endl; Print_Array_2D(lu[0], n, n);
#endif

	// Calc. x'

	double xp[n];
	xp[0]=b[0]/lu[0][0];
	for(int i=1; i<n; i++){
		double sum_lxp = 0.0;
		for(int j=0; j<i; j++) sum_lxp += lu[i][j]*xp[j];
		xp[i] = (b[i]-sum_lxp)/lu[i][i];
	}
#if _DEBUG_
	cout << "Vector x' is" << endl; Print_Array_2D(xp, n, 1);
#endif

	// Calc. x
	x[n-1]=xp[n-1];
	for(int i=n-2; i>=0; i--){
		double sum_ux = 0.0;
		for(int j=i+1; j<n; j++) sum_ux += lu[i][j]*x[j];
		x[i] = xp[i]-sum_ux;
	}

	return 1;
}


void Test_LU_Decomposition(){

/*
	const int n = 4;
	double a[n][n]={
			{ 0.0 , 1.0 , 3.0 ,-5.0 },
			{ 3.0 ,-2.0 , 1.0 , 1.0 },
			{ 8.0 ,-5.0 ,-7.0 , 2.0 },
			{ 1.0 , 1.0 ,-2.0 ,-3.0 }
	};
	double b[n]={-9.0, 6.0,	-15.0, -15.0};
	cout << "Matrix A is" << endl; Print_Array_2D(a[0], n, n);
	cout << "Vector b is" << endl; Print_Array_2D(b, n, 1);
	Pivot(a[0],b,n);
*/

// Test by Wolfram Alpha : x+2y+3z=14,-x+z=2,3x+3y+4z=21
	const int n = 3;
	double a[n][n]={
			{ 1.0, 2.0 , 3.0 },
			{ -1.0, 0.0, 1.0 },
			{ 3.0, 3.0, 4.0 }
	};
	double b[n]={14.0, 2.0, 21.0};

	double x[n];

	cout << "Matrix A is" << endl; Print_Array_2D(a[0], n, n);
	cout << "Vector b is" << endl; Print_Array_2D(b, n, 1);

	if(LU_Decomposition(a[0], n, b, x)>0){
		cout << "Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
	}else{
		cout << "Error" << endl;
	}

	if(LU_Decomposition2(a[0], n, b, x)>0){
		cout << "Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
	}else{
		cout << "Error" << endl;
	}

}

/*
 * 反復解法による連立１次方程式の解： ヤコビ法
 * 対角要素の行列A1と非対角要素の行列A2に分離して、第p近似解x(p)を
 *	x(p)=A1^-1*b-A1^-1*A2*x(0)
 * として求める。
*/
int JACOBI(double *a, double *b, double *x, int n, double eps, int kpp){

	double xx[n];

	///// Initialize
	for(int i=0; i<n; i++){
		xx[i] = 0.0;
		x[i] = 0.0;
		// Calc. Q
		b[i] /= a[i*n+i];
	}
	double t;
	for(int i=0; i<n; i++){
		t=a[i*n+i];
		for(int j=0; j<n; j++){
			// Calc. R
			a[i*n+j] /= t;
		}
	}

	///// Iteration for Solve
	int kp=0;
	int flag;
	do{
		kp++;
		flag=0;
		for(int i=0; i<n; i++){
			double sum = 0.0;
			for(int j=0; j<n; j++){
				// Calc. 2nd term of Eq.(2-4) or (2-10)
				if(i!=j) sum += xx[j]*a[i*n+j];
			}
			// Calc. 1st term of Eq.(2-4) or (2-10)
			x[i] = b[i]-sum;
		}
#if _DEBUG_
		// cout << "@kp=" << kp << ": Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
#endif
		/// Converse Judgement
		for(int i=0; i<n; i++) if(abs(x[i]-xx[i])>eps) flag = 1;

		// Prepare for next iteration
		for(int i=0; i<n; i++) xx[i] = x[i];
	}while(kp<kpp && flag>0);

	return kp;
}

/*
 * 反復解法による連立１次方程式の解： ガウス・ザイデル法
 * ヤコビ法の改良手法
 * ヤコビ法では反復pにおいて前の反復p-1での近似解xi(p-1)(i=1-n)を用いるが、
 * ガウス・ザイデル法では反復pにおいて既に計算済みのm(<n)個までのxiはxi(p-1)ではなくxi(p)(i=1-m)を用いる。
*/
int GaussSeidel(double *a, double *b, double *x, int n, double eps, int kpp){

	double xx[n];

	///// Initialize
	for(int i=0; i<n; i++){
		xx[i] = 0.0;
		x[i] = 0.0;
		// Calc. Q
		b[i] /= a[i*n+i];
	}
	double t;
	for(int i=0; i<n; i++){
		t=a[i*n+i];
		for(int j=0; j<n; j++){
			// Calc. R
			a[i*n+j] /= t;
		}
	}

	///// Iteration for Solve
	int kp=0;
	int flag;
	do{
		kp++;
		flag=0;
		for(int i=0; i<n; i++){
			double sum = 0.0;
			for(int j=0; j<n; j++){
				// Calc. 2nd term of Eq.(2-4) or (2-10)
				if(i!=j) sum += x[j]*a[i*n+j];
			}
			// Calc. 1st term of Eq.(2-4) or (2-10)
			x[i] = b[i]-sum;
		}
#if _DEBUG_
		 cout << "@kp=" << kp << ": Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
#endif
		/// Converse Judgement
		for(int i=0; i<n; i++) if(abs(x[i]-xx[i])>eps) flag = 1;

		// Prepare for next iteration
		for(int i=0; i<n; i++) xx[i] = x[i];
	}while(kp<kpp && flag>0);

	return kp;
}


int SOR(double *a, double *b, double *x, int n, double eps, int kpp, double w=0.0){

	double xx[n];

	if(w==0.0){
		double Pai=3.141516;
		double s=Pai/n;
		w=2.0/(1.0+sin(s));
		cout << "Relaxation Parameter w = " << w << endl;
	}

	///// Initialize
	for(int i=0; i<n; i++){
		xx[i] = 0.0;
		x[i] = 0.0;
		// Calc. Q
		b[i] /= a[i*n+i];
	}
	double t;
	for(int i=0; i<n; i++){
		t=a[i*n+i];
		for(int j=0; j<n; j++){
			// Calc. R
			a[i*n+j] /= t;
		}
	}

	///// Iteration for Solve
	int kp=0;
	int flag;
	do{
		kp++;
		flag=0;
		for(int i=0; i<n; i++){
			double sum = 0.0;
			for(int j=0; j<n; j++){
				// Calc. 2nd term of Eq.(2-4) or (2-10)
				// if(i!=j) sum += x[j]*a[i*n+j];
				sum += x[j]*a[i*n+j];
			}
			// Calc. 1st term of Eq.(2-4) or (2-10)
			x[i] = b[i]-sum;
			x[i] = x[i]*w + xx[i];
		}
#if _DEBUG_
		 cout << "@kp=" << kp << ": Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
#endif
		/// Converse Judgement
		for(int i=0; i<n; i++) if(abs(x[i]-xx[i])>eps) flag = 1;

		// Prepare for next iteration
		for(int i=0; i<n; i++) xx[i] = x[i];
	}while(kp<kpp && flag>0);

	return kp;
}


void Test_JACOBI_GaussSeidel_SOR(){

	int kp, kpp;

#define _TEST_SET_NUM_SOR_ 7
#if _TEST_SET_NUM_SOR_ == 1
	// Test by Wolfram Alpha : 7w+2x+-1y+1z=12,1w+5x+1y+-2z=6,2w+3x+8y+1z=36,2w+-2x+-1y+10z=35
	const int n = 4;
	double a[n][n]={
			{7.0, 2.0,-1.0, 1.0},
			{1.0, 5.0, 1.0,-2.0},
			{2.0, 3.0, 8.0, 1.0},
			{2.0,-2.0,-1.0,10.0}
	};
	double b[n]={12.0, 6.0,36.0,35.0};
	kpp=50;
#elif _TEST_SET_NUM_SOR_ == 2
	const int n = 3;
	double a[n][n]={
			{10.0, 1.0, 1.0},
			{ 2.0,10.0, 2.0},
			{ 1.0, 2.0,10.0}
	};
	double b[n]={15.0, 28.0,35.0};
#elif _TEST_SET_NUM_SOR_ == 3
	const int n = 3;
	double a[n][n]={
			{ 1.0, 2.0,-2.0},
			{ 1.0, 1.0, 1.0},
			{ 2.0, 2.0, 1.0}
	};
	double b[n]={ 2.0, 4.0, 6.0};
#elif _TEST_SET_NUM_SOR_ == 4
	const int n = 2;
	double a[n][n]={
			{ 1.0, 3.0 },
			{ 2.0, 1.0 }
	};
	double b[n]={ 7.0, 3.0 };
#elif _TEST_SET_NUM_SOR_ == 5
	const int n = 2;
	double a[n][n]={
			{ 5.0, 1.0 },
			{ 2.0, 5.0 }
	};
	double b[n]={ 6.0, 7.0 };
#elif _TEST_SET_NUM_SOR_ == 6
	const int n = 3;
	double a[n][n]={
			{10.0, 3.0, 1.0},
			{-2.0, 5.0, 1.0},
			{ 2.0, 1.0, 5.0}
	};
	double b[n]={ 14.0, 4.0, 8.0};
#elif _TEST_SET_NUM_SOR_ == 7
	const int n = 2;
	double a[n][n]={
			{ 1.0, 2.0 },
			{ 3.0, 8.0 }
	};
	double b[n]={ 4.0, 14.0 };
#endif


	cout << "Input Maximum Iteration (kpp) : " << endl; cin >> kpp;


	double x[n];

	cout << "Matrix A is" << endl; Print_Array_2D(a[0], n, n);
	cout << "Check Diagonal Priority, Mi(<1.0) is " << Check_DiagonalPriority(a[0],n) << endl;
	cout << "Vector b is" << endl; Print_Array_2D(b, n, 1);


	cout << "By JACOBI" << endl;
	kp = JACOBI(a[0], b, x, n, 0.000001, kpp);
	if(kp<kpp){
		cout << "Converged@kp=" << kp << endl;
		cout << "Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
	}else{
		cout << "Error" << endl;
	}

	cout << "By GaussSeidel" << endl;
	kp = GaussSeidel(a[0], b, x, n, 0.000001, kpp);
	if(kp<kpp){
		cout << "Converged@kp=" << kp << endl;
		cout << "Solution Vector x is" << endl; Print_Array_2D(x, n, 1);
	}else{
		cout << "Error" << endl;
	}

	cout << "By SOR" << endl;
	kp = SOR(a[0], b, x, n, 0.000001, kpp);
	if(kp<kpp){
		cout << "Converged@kp=" << kp << endl;
		cout << "Solution Vectorx is" << endl; Print_Array_2D(x, n, 1);
	}else{
		cout << "Error" << endl;
	}

}


double Interpolate_Lagrange_Base(double xn[], double fn[], int n, double x){

	double Qi;
	double Pn=0.0;
	for(int i=0; i<n ; i++){
		double Qi_Numerator = 1.0;
		double Qi_Denominator=1.0;
		for(int j=0; j<n; j++){
			if (i!=j){
				Qi_Numerator *= (x-xn[j]);
				Qi_Denominator *= (xn[i]-xn[j]);
			}
		}

		Qi=Qi_Numerator/Qi_Denominator;
#if _DEBUG_
			cout << "fn[" << i << "]=" << fn[i] <<", Qi=" << Qi << endl;
#endif
		Pn += fn[i]*Qi;
	}

	return(Pn);
}


double Interpolate_Lagrange(double xn[], double fn[], int n, double x){

	int k=-1;
	double Qi_Numerator0=1.0;
	for(int i=0; i<n; i++){
		if(x==xn[i]){
			k=i;
#if _DEBUG_
			cout << "k=" << k << endl;
#endif
		}
		Qi_Numerator0 *= (x-xn[i]);
	}

	double Qi;
	double Pn=0.0;
	for(int i=0; i<n ; i++){
		if(k==-1){
			double Qi_Numerator = Qi_Numerator0/(x-xn[i]);

			double Qi_Denominator=1.0;
			for(int j=0; j<n; j++){
				if (i!=j){
					Qi_Denominator *= (xn[i]-xn[j]);
				}
			}

			Qi=Qi_Numerator/Qi_Denominator;
		}else{
			if(k==i){
				Qi=1.0;
			}else{
				Qi=0.0;
			}
		}
#if _DEBUG_
			cout << "fn[" << i << "]=" << fn[i] <<", Qi=" << Qi << endl;
#endif
		Pn += fn[i]*Qi;
	}

	return(Pn);
}


void Test_Interpolate_Lagrange(){

	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};

	double x;
	cout << "Input x" << endl; cin >> x;

	double Pn;
	cout << "Interpolate by Lagrange (Basic Algorithm)" << endl;
	Pn=Interpolate_Lagrange_Base(xn, fn, n, x);
	cout << "Pn(" << x << ")= " << Pn << endl;

	cout << "Interpolate by Lagrange (Light Algorithm)" << endl;
	Pn=Interpolate_Lagrange(xn, fn, n, x);
	cout << "Pn(" << x << ")= " << Pn << endl;

}


double GetValue_Polynominal(double a[], int n, double x){

	double f=0.0;
	for(int i=0; i<n; i++){
		f += a[i]*pow(x,i);
	}
	return f;
}

//
// (an*x^n + ... + a2*x^2 + a1*x^1 + a0*x^0)*(x+b)
double *Get_Polynominal(double *a, int n, double b){

#if _DEBUG_ && _DEBUG_ < 2
	Print_Array(a,n);
#endif

	a = (double *)realloc(a,sizeof(double) * (n+1));
//	Print_Array_2D(a,1,n+1);

	a[n] = a[n-1];

//	Print_Array(a,n+1);

	for(int i=n-1; i>0; i--){
		// a[i] = a[i-1]+a[i]*b;
		a[i]*=b;
		a[i]+=a[i-1];
//		Print_Array(a,n+1);
	}
	a[0]*=b;
#if _DEBUG_ && _DEBUG_ < 2
	Print_Array(a,n+1);
#endif

	return a;
}


void Test_Get_Polynominal(){

	int n;
	cout << "Input n" << endl; cin >> n;

	double *a = (double *)malloc(sizeof(double) * n);

	for(int i=0; i<n; i++){
		cout << "Input a[" << i << " ]=" << endl; cin >> a[i];
	}
//	Print_Array(a,n);

	double b;
	cout << "Input b" << endl; cin >> b;
	a=Get_Polynominal(a,n,b);
	n++;
	Print_Array(a,n);

	cout << "Input b" << endl; cin >> b;
	a=Get_Polynominal(a,n,b);
	n++;
	Print_Array(a,n);

	free(a);

	return;
}


void Interpolate_Lagrange(double xn[], double fn[], int n, double Pn[]){



	for (int i = 0; i < n; ++i) {
		int nn=1;
		double *Qi_Numerator = (double *)malloc(sizeof(double) * nn);
		Qi_Numerator[0]=1.0;
		double Qi_Denominator = 1.0;
		for (int j = 0; j < n; ++j) {
			if (i!=j) {
				Qi_Numerator=Get_Polynominal(Qi_Numerator, nn++, -1.0*xn[j]);
				Qi_Denominator *= (xn[i]-xn[j]);
			}
		#if _DEBUG_ && _DEBUG_ < 3
			cout << "Qi_Numerator:"; Print_Array(Qi_Numerator,nn); cout << endl;
			cout << "Qi_Denominator=" << Qi_Denominator << endl;
		#endif
		}
	#if _DEBUG_ && _DEBUG_ < 4
		cout << "Qi_Numerator:"; Print_Array(Qi_Numerator,nn);
		cout << "Qi_Denominator=" << Qi_Denominator << endl << endl;
	#endif
		free(Qi_Numerator);
	}

	return;
}

void Test_Interpolate_Lagrange2(){

	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	double Pn[] = {0.0, 0.0, 0.0, 0.0};



	cout << "Interpolate by Lagrange (Light Algorithm)" << endl;
	Interpolate_Lagrange(xn, fn, n, Pn);
//	Print_Array(Pn,n);

	return;
}

double Function_for_GetDerivertiveValueBy5PointsFormula(double x){
	return exp(x);
}

double GetDerivertiveValueBy5PointsFormula(double (*func)(double), double a, double h){
	return (func(a-2.0*h) - 8.0*func(a-1.0*h) + 8.0*func(a+h) - func(a+2.0*h))/12.0/h;
}

void Test_GetDerivertiveValueBy5PointsFormula(){
	double val1stDeriv = GetDerivertiveValueBy5PointsFormula(Function_for_GetDerivertiveValueBy5PointsFormula,1,0.1);
	cout << "Derivertive of Func. = " << val1stDeriv << endl;
}

double Function_for_GetIntegralValue(double x){
//	return sqrt(1.0-x*x);
	return 1.0/(1.0+x*x);
//	return x*x*log(x);
}

double GetIntegralValueByTrapezoidalMethod(double (*func)(double),	double x1 ,double x2,int numDiv){

	double x, f;
	double h = (x2 - x1)/numDiv;
	double sum = 0.0;

	for (int i = 0; i <= numDiv; ++i) {
		x = x1+i*h;
		f =func(x);
//		cout << "Value of Function at " << x << " = " << f << endl;
		sum += ((i==0)||(i==numDiv)? 1.0 : 2.0) * f;
	}

	return h/2.0*sum;
}

// incr : >0 = numDiv += incr, <0 = numDiv *= -incr
double GetIntegralValueByIterativeTrapezoidalMethod(double (*func)(double), double x1 ,double x2,int numDiv, double RelErr, int incr=10){

	double Sigma0, Sigma;
	int i=1;
	Sigma= GetIntegralValueByTrapezoidalMethod(func, x1, x2, numDiv);
	cout << i << " : Integral Value = " << Sigma << endl;
	Sigma0 = Sigma;
	do {
		++i;
		if(incr<0){
			numDiv *= -incr;
		}else{
			numDiv += incr;
		}
		Sigma0 = Sigma;
		Sigma= GetIntegralValueByTrapezoidalMethod(func, x1, x2, numDiv);
		cout << i << " : Integral Value = " << Sigma << endl;

		RelErr = abs((Sigma - Sigma0)/Sigma0);
		cout << i << " : Division Number = " << numDiv << ", Relative Error = " << RelErr << endl;
	} while (RelErr>1e-6);

	return Sigma;
}

double GetIntegralValueBySimpsonsFormula(double (*func)(double), double x1 ,double x2, int numDiv){

	double x, f;

	if((numDiv % 2) != 0) numDiv++;
	double h = (x2 - x1)/numDiv;

	double sum = 0.0;
	sum += func(x1);
	for (int i = 1; i < numDiv; ++i) {
		x = x1+i*h;
		f =func(x);
//		cout << "Value of Function at " << x << " = " << f << endl;
		if((i % 2) == 0){
			sum += 2.0*f;
		}else{
			sum += 4.0*f;
		}
	}
	sum += func(x2);

	return h/3.0*sum;
}
/*
double GetIntegralValueByLagrangeInterpolation(double (*func)(double), double x1 ,double x2, int numDiv){

	double h = (x2 - x1)/numDiv;
	int n = numDiv+1;
	double xn[n];
	double fn[n];

	for (int i = 0; i < n; ++i) {
		xn[i] = x1 + i*h;
		fn[i] = func(xn[i]);
	}
	Print_Array(xn,n);
	Print_Array(fn,n);

	Polynominal Pn = Polynominal::SetAsLagrangeInterpolation(xn,fn,n); Pn.Show('x',-1);
	Polynominal IntegPn = Pn.GetIntegral(); IntegPn.Show('x',-1);

	return IntegPn.GetValue(x2) - IntegPn.GetValue(x1);

}
*/
/*
void Test_GetIntegralValue(){

	double Sigma;
	double x1 = 0.0;
	double x2 = 1.0;
	int numDiv = 10;

	cout << "Iterative Trapezoidal Method" << endl;
//	double Sigma = GetIntegralValueByIterativeTrapezoidalMethod(Function_for_GetIntegralValue,x1,x2,numDiv,1e-6, 10);
	Sigma = GetIntegralValueByIterativeTrapezoidalMethod(Function_for_GetIntegralValue,x1,x2,numDiv,1e-6, -2);
	cout << "Integral Value = " << Sigma << endl;
	cout << "Pi = " << Sigma*4.0 << endl;

	cout << "Simpsons Formula" << endl;
	Sigma = GetIntegralValueBySimpsonsFormula(Function_for_GetIntegralValue,x1,x2,numDiv);
	cout << "Integral Value = " << Sigma << endl;
	cout << "Pi = " << Sigma*4.0 << endl;

	cout << "LagrangeInterpolation" << endl;
	Sigma = GetIntegralValueByLagrangeInterpolation(Function_for_GetIntegralValue,x1,x2,numDiv);
	cout << "Integral Value = " << Sigma << endl;
	cout << "Pi = " << Sigma*4.0 << endl;
}
*/



int main() {
	cout << "!!! Welcome to Introduction and Applications to Numerical Calculations by C++!!!" << endl;

//	Test_GaussJordan();
//	Test_LU_Decomposition();
//	Test_InverseMatrixByGaussJordan();
	Test_JACOBI_GaussSeidel_SOR();
//	Test_Interpolate_Lagrange();
//	Test_Get_Polynominal();
//	Test_Interpolate_Lagrange2();
//	Test_Class_DividedDifferencelTable();
//	Test_Class_Polynominal();
//	Test_GetDerivertiveValueBy5PointsFormula();
//	Test_GetIntegralValue();
//	Test_Class_GaussIntegralFormula();
	return 1;


}
