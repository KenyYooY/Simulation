/*
 * Polynominal.cpp
 *
 *  Created on: 2016/05/26
 *      Author: Kanie
 */

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

class DividedDifferencelTable{
private:
	int n;	// Number of Points ((x0,f0), (x1,f1), ... , (xn-1,fn-1)
	double *xn;
	int nt; // Number of Table
	double *fn;

public:
	DividedDifferencelTable();
	DividedDifferencelTable(double *, int);
	DividedDifferencelTable(double *, double *, int);
	~DividedDifferencelTable(){
//		free(this->xn);
		if( this->xn != NULL) { delete[] this->xn; }
//		free(this->fn);
		if( this->fn != NULL ){ delete[] this->fn; }
		this->n = 0;
		this->nt = 0;
	};
	int GetNumberOfPoints(){ return this->n; };
	int GetSizeOfTable(){return this->nt; };
	double GetTableValue(int, int);
	double GetTableValueForward(int);
	double GetTableValueBackward(int);
	double GetTableValueCentralForward(int);
	double GetTableValueCentralBackward(int);
	double GetTableValue(int);
	void ShowTable(int, char);
	void ShowAllTable(char);

};

// dim = 0, ... , n-1
void DividedDifferencelTable::ShowTable(int dim, char delimiter=','){
	if( (dim >= 0) && (dim < this->n) ){
		int start = 0;
		for (int i = 0; i < this->n; ++i) {
			int j_max = this->n - i;
			if(dim == i){
				for (int j = 0; j < j_max; ++j) {
					cout << this->fn[start+j] << delimiter;
				}
				cout << endl;
				return;
			}
			start += j_max;
		}
	}else{
		cout << "Error : Out of Range in Input \"dim\"";
	}

	return;
}

void DividedDifferencelTable::ShowAllTable(char delimiter=','){

	int start = 0;
	for (int i = 0; i < this->n; ++i) {
		int j_max = this->n - i;
		for (int j = 0; j < j_max; ++j) {
			cout << this->fn[start+j] << delimiter;
		}
		cout << endl;
		start += j_max;
	}

	return;
}

DividedDifferencelTable::DividedDifferencelTable(){
	this->n = 0;
	this->xn = NULL;
	this->nt = 0;
	this->fn = NULL;
}

#define _DEBUG_DividedDifferencelTable_DividedDifferencelTable_1_ 0
DividedDifferencelTable::DividedDifferencelTable(double fn_set[], int n_set){

	this->n = n_set;
	this->xn = NULL;
	this->nt = 0;
	for (int i = 0; i < this->n; ++i) {
		this->nt += i+1;
	}
//	this->fn = (double *)malloc(sizeof(double) * this->nt);
	this->fn = new double[this->nt];


	// Set 0th Order
	for (int i = 0; i < this->n; ++i) {
		this->fn[i] = fn_set[i];
	}

	// Set 1st~ Order
	int start_pre = 0;
	int start = this->n;
	for (int i = 1; i < this->n; ++i) {
		int j_max = this->n - i;
		for (int j = 0; j < j_max; ++j) {
			this->fn[start+j] = ( this->fn[start_pre+j+1] - this->fn[start_pre+j] );
		#if _DEBUG_DividedDifferencelTable_DividedDifferencelTable_1_ > 1
				cout << "start+j = " << start+j << endl;
				cout << "start_pre+j+1 = " << start_pre+j+1 << ", start_pre+j = " << start_pre+j << endl;
				cout << "j+i = " << j+i << ", j = " << j << endl;
		#endif
		}
		start_pre = start;
		start += j_max;
	#if _DEBUG_DividedDifferencelTable_DividedDifferencelTable_1_ > 0
		cout << "start_pre = " << start_pre << endl;
		cout << "start = " << start << endl;
	#endif
	}
}


#define _DEBUG_DividedDifferencelTable_DividedDifferencelTable_2_ 0
DividedDifferencelTable::DividedDifferencelTable(double xn_set[], double fn_set[], int n_set){

	this->n = n_set;
//	this->xn = (double *)malloc(sizeof(double) * this->n);
	this->xn = new double[this->n];
	this->nt = 0;
	for (int i = 0; i < this->n; ++i) {
		this->nt += i+1;
	}
//	this->fn = (double *)malloc(sizeof(double) * this->nt);
	this->fn = new double[this->nt];

	// Set 0th Order
	for (int i = 0; i < this->n; ++i) {
		this->xn[i] = xn_set[i];
		this->fn[i] = fn_set[i];
	}

	// Set 1st~ Order
	int start_pre = 0;
	int start = this->n;
	for (int i = 1; i < this->n; ++i) {
		int j_max = this->n - i;
		for (int j = 0; j < j_max; ++j) {
			this->fn[start+j] = ( this->fn[start_pre+j+1] - this->fn[start_pre+j] ) / ( this->xn[j+i] - this->xn[j]);
		#if _DEBUG_DividedDifferencelTable_DividedDifferencelTable_2_ > 1
				cout << "start+j = " << start+j << endl;
				cout << "start_pre+j+1 = " << start_pre+j+1 << ", start_pre+j = " << start_pre+j << endl;
				cout << "j+i = " << j+i << ", j = " << j << endl;
		#endif
		}
		start_pre = start;
		start += j_max;
	#if _DEBUG_DividedDifferencelTable_DividedDifferencelTable_2_ > 0
		cout << "start_pre = " << start_pre << endl;
		cout << "start = " << start << endl;
	#endif
	}
}

// order = 0, ... , n-1, index = 0, ... , n-1
double DividedDifferencelTable::GetTableValue(int order, int index){
	if( (order >= 0) && (order < this->n) ){
		int start = 0;
		for (int i = 0; i < this->n; ++i) {
			int j_max = this->n - i;
			if(order == i){
				if((index>=0) && (index<j_max)){
					for (int j = 0; j < j_max; ++j) {
						if(index == j){
							return this->fn[start+j];
						}
					}
				}else{
					cout << "Error : Out of Range in Input \"index\" :" << index << endl;
					break;
				}
			}
			start += j_max;
		}
	}else{
		cout << "Error : Out of Range in Input \"order\" :" << index << endl;
	}
	return 0.0;
}

// order = 0, ... , n-1
double DividedDifferencelTable::GetTableValueForward(int order){
	return this->GetTableValue(order,0);
}

// order = 0, ... , n-1
double DividedDifferencelTable::GetTableValueBackward(int order){
	return this->GetTableValue(order,this->n - order - 1);
}

// order = 0, ... , n-1
double DividedDifferencelTable::GetTableValueCentralForward(int order){

	int index = 0;
	int nTableAtOrder = this->n - order;

	if((this->n % 2) == 0){
		if((nTableAtOrder % 2) == 0){
			index = (int)(nTableAtOrder/2)-1;
		}else{
			index = (int)(nTableAtOrder/2);
		}
	}else{
		index = (int)(nTableAtOrder/2);
	}

	return this->GetTableValue(order,index);

}

// order = 0, ... , n-1
double DividedDifferencelTable::GetTableValueCentralBackward(int order){

	int index = 0;
	int nTableAtOrder = this->n - order;

	if((this->n % 2) == 0){
		index = (int)(nTableAtOrder/2);
	}else{
		if((nTableAtOrder % 2) == 0){
			index = (int)(nTableAtOrder/2) - 1;
		}else{
			index = (int)(nTableAtOrder/2);
		}
	}

	return this->GetTableValue(order,index);

}


void Test_Class_DividedDifferencelTable(){

#define _DEBUG_TEST_CLASS_DividedDifferencelTable_ 2

#if _DEBUG_TEST_CLASS_DividedDifferencelTable_ == 1
	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	DividedDifferencelTable dffTab = DividedDifferencelTable(xn,fn,n);

	int nn = dffTab.GetNumberOfPoints();
	cout << "Number of Points = " << nn << ", Size of Table = " << dffTab.GetSizeOfTable() << endl;
	dffTab.ShowAllTable();
	for (int i = 0; i < nn; ++i) {
		dffTab.ShowTable(i);
	}
	int start = 0;
	for (int i = 0; i < nn; ++i) {
		int j_max = nn - i;
		for (int j = 0; j < j_max; ++j) {
			cout << dffTab.GetTableValue(i,j) << ",";
		}
		cout << endl;
		start += j_max;
	}
#elif _DEBUG_TEST_CLASS_DividedDifferencelTable_ == 2
	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	DividedDifferencelTable dffTab = DividedDifferencelTable(xn,fn,n);
	dffTab.ShowAllTable();
	int nn = dffTab.GetNumberOfPoints();

	cout << "Table Values for Newton Forward Interpolation" << endl;
	for (int i = 0; i < nn; ++i) {
		double tmp = dffTab.GetTableValueForward(i);
		cout << tmp << ", ";
	}
	cout << endl;

	cout << "Table Values for Newton Backward Interpolation" << endl;
	for (int i = 0; i < nn; ++i) {
		double tmp = dffTab.GetTableValueBackward(i);
		cout << tmp << ", ";
	}
	cout << endl;

	cout << "Table Values for Gauss Forward Interpolation" << endl;
	for (int i = 0; i < nn; ++i) {
		double tmp = dffTab.GetTableValueCentralForward(i);
		cout << tmp << ", ";
	}
	cout << endl;

	cout << "Table Values for Gauss Backward Interpolation" << endl;
	for (int i = 0; i < nn; ++i) {
		double tmp = dffTab.GetTableValueCentralBackward(i);
		cout << tmp << ", ";
	}
	cout << endl;

#endif

	return;
}


// f(x) = a0 + a1*x + a2*x^2 + ... + an-1*x^(n-1)
class Polynominal{
private:
//	int basis;
	int n;
	double *an;

	static int Factorial(int);

public:
	Polynominal();
	Polynominal(const Polynominal&); // Copy Constructor for Initialize
	Polynominal(int);
	Polynominal(double);
	Polynominal(double, double);
	Polynominal(int, double *);
	static Polynominal SetAsLagrangeInterpolation(double [], double [], int);
	static Polynominal SetAsHermiteInterpolation(double [], double [], double [], int);
	static Polynominal SetAsNewtonInterpolation(double [], double [], int);
	static Polynominal SetAsNewtonForwardDifferenceInterpolation(double [], int);
	static Polynominal SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation(double [], int, double);
	static Polynominal SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation(double [], int, double);
	static Polynominal SetAsNewtonBackwardDifferenceInterpolation(double [], int);
	static Polynominal SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation(double [], int, double);
	static Polynominal SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation(double [], int, double);
	static Polynominal SetAsGaussForwardInterpolation(double [], int);
	static Polynominal SetAs1stOrderDerivertiveByGaussForwardInterpolation(double [], int, double);
	static Polynominal SetAs2ndOrderDerivertiveByGaussForwardInterpolation(double [], int, double);
	static Polynominal SetAsGaussBackwardInterpolation(double [], int);
	static Polynominal SetAs1stOrderDerivertiveByGaussBackwardInterpolation(double [], int, double);
	static Polynominal SetAs2ndOrderDerivertiveByGaussBackwardInterpolation(double [], int, double);
	static Polynominal SetAsStirlingInterpolation(double [], int);
	static Polynominal SetAsLegendre(int);
	static Polynominal SetAsHermite(int);
	static Polynominal SetAsLaguerre(int);
	~Polynominal(){
//		free(this->an);
		if(this->an != NULL){ delete[] this->an; }
		this->n=0;
	};
	Polynominal CheckMaxOrder();
	Polynominal& operator=(const Polynominal&);
	Polynominal operator+(Polynominal);
	friend Polynominal operator+(Polynominal, double);
	friend Polynominal operator+(double, Polynominal);
	Polynominal operator-(Polynominal);
	friend Polynominal operator-(Polynominal, double);
	friend Polynominal operator-(double, Polynominal);
	Polynominal operator*(Polynominal);
	friend Polynominal operator*(Polynominal, double);
	friend Polynominal operator*(double, Polynominal);
//	friend Polynominal SyntheticDivision(Polynominal, Polynominal, Polynominal *);
	friend Polynominal operator/(Polynominal, double);
	Polynominal Power(unsigned int);
	Polynominal Substitute(Polynominal);
	Polynominal GetDerivertive();
	Polynominal GetIntegral();
	void Show(char, int);
	int GetDimenstion(){ return this->n; }
	double GetValue(double);
};

Polynominal Polynominal::CheckMaxOrder(){

	int MaxOrder;

	for (int i = 0; i < this->n; i++) {
		if(this->an[this->n-1-i]!=0.0){
			MaxOrder = this->n-i;
			break;
		}
	}

	Polynominal ret(MaxOrder);
	for (int i = 0; i < MaxOrder; ++i) {
		ret.an[i]=this->an[i];
	}

	return ret;
}

void Polynominal::Show(char delimiter = ',', int direc = 1){

	int ii;

	for(int i=0; i<this->n; i++){
		ii = direc>=0 ? i : this->n-i-1;
		if(isalpha(delimiter)){
			cout << this->an[ii] << '*' << delimiter << '^' << ii;
			if(i<this->n-1){ cout << " + "; }
		}else{
			cout << this->an[ii];
			if(i<this->n-1){ cout << delimiter; }
		}
	}
	cout << endl;

}


Polynominal::Polynominal(){
	this->n=0;
	this->an = NULL;
}

// Copy Constructor for Initialize
Polynominal::Polynominal(const Polynominal &Pn){
	this->n = Pn.n;
//	this->an = (double *)malloc(sizeof(double) * this->n);
	this->an = new double[this->n];
	for (int i = 0; i < this->n; ++i) {
		this->an[i] = Pn.an[i];
	}
}

// f(x) = a0 + a1*x + a2*x^2 + ... + an-1*x^(n-1), ai=0.0(i=0...n-1)
Polynominal::Polynominal(int n_set){
//	if(n>0) free(an);
	this->n=n_set;
//	this->an = (double *)malloc(sizeof(double) * this->n);
	this->an = new double[this->n];

	for (int i = 0; i < this->n; ++i) {
		this->an[i]=0.0;
	}

}

// f(x) = a0
Polynominal::Polynominal(double a0){
//	if(n>0) free(an);
	this->n = 1;
//	this->an = (double *)malloc(sizeof(double) * this->n);
	this->an = new double[this->n];
	an[0]=a0;
}

// f(x) = a0 + a1*x
Polynominal::Polynominal(double a0, double a1){
//	if(n>0) free(an);
	this->n=2;
//	this->an = (double *)malloc(sizeof(double) * this->n);
	this->an = new double[this->n];
	this->an[0]=a0;
	this->an[1]=a1;
}

// f(x) = a0 + a1*x + a2*x^2 + ... + an-1*x^(n-1)
Polynominal::Polynominal(int n_set, double an_set[]){
//	if(n>0) free(an);
	this->n = n_set;
//	this->an = (double *)malloc(sizeof(double) * this->n);
	this->an = new double[this->n];

	for (int i = 0; i < this->n; ++i) {
		this->an[i] = an_set[i];
	}

}

#define _DEBUG_Polynominal_SetAsLagrangeInterpolation_ 0
Polynominal Polynominal::SetAsLagrangeInterpolation(double xn[], double fn[], int n_set){

	Polynominal Pn(0.0);

	for (int i = 0; i < n_set; ++i) {
		Polynominal Qi(1.0);
		double Qi_Denominator = 1.0;
		for (int j = 0; j < n_set; ++j) {
			if (i!=j) {
				Qi = Qi * Polynominal(-1.0*xn[j],1);
				Qi_Denominator *= (xn[i]-xn[j]);
			}
		#if _DEBUG_Polynominal_SetAsLagrangeInterpolation_ > 1
			cout << "Qi:"; Qi.Show();
			cout << "Qi_Denominator=" << Qi_Denominator << endl;
		#endif
		}
		Qi = Qi/Qi_Denominator;
		Pn = Pn + fn[i]*Qi;
	#if _DEBUG_Polynominal_SetAsLagrangeInterpolation_ > 0
		cout << "Qi:"; Qi.Show('x',-1);
		Pn.Show('x',-1); cout << endl;
	#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsHermiteInterpolation_ 0
Polynominal Polynominal::SetAsHermiteInterpolation(double xn[], double fn[], double fpn[], int n_set){

	Polynominal Pn(0.0);

	for (int i = 0; i < n_set; ++i) {
		Polynominal Qi(1.0);
		Polynominal Qpi;
		double Qi_Denominator = 1.0;
		for (int j = 0; j < n_set; ++j) {
			if (i!=j) {
				Qi = Qi * Polynominal(-1.0*xn[j],1);
				Qi_Denominator *= (xn[i]-xn[j]);
			}
		#if _DEBUG_Polynominal_SetAsHermiteInterpolation_ > 1
			cout << "Qi_Numerator:"; Qi.Show();
			cout << "Qi_Denominator=" << Qi_Denominator << endl;
		#endif
		}
		Qi = Qi / Qi_Denominator;
		Qpi = Qi.GetDerivertive();
		Polynominal Qi2 = Qi*Qi;
		Polynominal Gi = (1.0 - 2.0*Polynominal(-1.0*xn[i],1.0)*Qpi)*Qi2;
		Polynominal Ri = Polynominal(-1.0*xn[i],1.0)*Qi2;
		Pn = Pn + fn[i]*Gi + fpn[i]*Ri;
	#if _DEBUG_Polynominal_SetAsHermiteInterpolation_ > 0
		cout << "Qi:"; Qi.Show('x',-1);
		cout << "Qi^2:"; Qi2.Show('x',-1);
		cout << "Qpi:"; Qpi.Show('x',-1);
		cout << "Gi:"; Gi.Show('x',-1);
		cout << "Ri:"; Ri.Show('x',-1);
		Pn.Show(); cout << endl;
	#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsNewtonInterpolation_ 0
Polynominal Polynominal::SetAsNewtonInterpolation(double xn[], double fn[], int n_set){

	Polynominal Pn(0.0);

	DividedDifferencelTable dffTab = DividedDifferencelTable(xn, fn, n_set);
#if _DEBUG_Polynominal_SetAsNewtonInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal tmp(1.0);
	for (int i = 0; i < n_set; ++i) {
		Pn = Pn + tmp*dffTab.GetTableValue(i,0);

		tmp = tmp*Polynominal(-1.0*xn[i],1.0);
	#if _DEBUG_Polynominal_SetAsNewtonInterpolation_ > 0
		tmp.Show('x',-1);
	#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsNewtonForwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAsNewtonForwardDifferenceInterpolation(double fn[], int n_set){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAsNewtonForwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

//	Polynominal Pn(dffTab.GetTableValue(0,0));
	Polynominal Pn(dffTab.GetTableValueForward(0));
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(1.0-i),1.0);
#if _DEBUG_Polynominal_SetAsNewtonForwardDifferenceInterpolation_ > 1
		tmp.Show('x',-1);
//		cout << "Table Value = " << dffTab.GetTableValue(i,0) << ", Factorial = " << Polynominal::Factorial(i) << endl;
		cout << "Table Value = " << dffTab.GetTableValueForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
//		Pn = Pn + tmp*dffTab.GetTableValue(i,0)/Polynominal::Factorial(i);
		Pn = Pn + tmp*dffTab.GetTableValueForward(i)/Polynominal::Factorial(i);

#if _DEBUG_Polynominal_SetAsNewtonForwardDifferenceInterpolation_ > 0
		Pn.Show();
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueForward(0));
	Pn = Pn.GetDerivertive();
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(1.0-i),1.0);
		Polynominal tmp1 = tmp.GetDerivertive();
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 1
		cout << "tmp : "; tmp.Show('x',-1);
		cout << "tmp1 : "; tmp1.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + 1.0/h*tmp1/Polynominal::Factorial(i)*dffTab.GetTableValueForward(i);

#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 0
		Pn.Show();
#endif
	}

	return Pn;

}

#define _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueForward(0));
	Pn = Pn.GetDerivertive();
	Pn = Pn.GetDerivertive();
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(1.0-i),1.0);
		Polynominal tmp1 = tmp.GetDerivertive();
		Polynominal tmp2 = tmp1.GetDerivertive();
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 1
		cout << "tmp : "; tmp.Show('x',-1);
		cout << "tmp1 : "; tmp1.Show('x',-1);
		cout << "tmp2 : "; tmp2.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + 1.0/h/h*tmp2/Polynominal::Factorial(i)*dffTab.GetTableValueForward(i);

#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation_ > 0
		Pn.Show();
#endif
	}

	return Pn;

}

#define _DEBUG_Polynominal_SetAsNewtonBackwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAsNewtonBackwardDifferenceInterpolation(double fn[], int n_set){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAsNewtonBackwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

//	Polynominal Pn(dffTab.GetTableValue(0,n_set-1));
	Polynominal Pn(dffTab.GetTableValueBackward(0));
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(i-1.0),1.0);
#if _DEBUG_Polynominal_SetAsNewtonBackwardDifferenceInterpolation_ > 1
		tmp.Show('x',-1);
//		cout << "Table Value = " << dffTab.GetTableValue(i,n_set-i-1) << ", Factorial = " << Polynominal::Factorial(i) << endl;
		cout << "Table Value = " << dffTab.GetTableValueBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
//		Pn = Pn + tmp*dffTab.GetTableValue(i,n_set-i-1)/Polynominal::Factorial(i);
		Pn = Pn + tmp*dffTab.GetTableValueBackward(i)/Polynominal::Factorial(i);

#if _DEBUG_Polynominal_SetAsNewtonBackwardDifferenceInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueBackward(0));
	Pn = Pn.GetDerivertive();
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(i-1.0),1.0);
		Polynominal tmp1 = tmp.GetDerivertive();
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 1
		cout << "tmp : " << tmp.Show('x',-1);
		cout << "tmp1 : " << tmp1.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h)*tmp1/Polynominal::Factorial(i)*dffTab.GetTableValueBackward(i);

#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ 0
Polynominal Polynominal::SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueBackward(0));
	Pn = Pn.GetDerivertive();
	Pn = Pn.GetDerivertive();
	Polynominal tmp(1.0);
	for (int i = 1; i < n_set; ++i) {
		tmp = tmp*Polynominal((double)(i-1.0),1.0);
		Polynominal tmp1 = tmp.GetDerivertive();
		Polynominal tmp2 = tmp1.GetDerivertive();
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 1
		cout << "tmp : " << tmp.Show('x',-1);
		cout << "tmp1 : " << tmp1.Show('x',-1);
		cout << "tmp2 : " << tmp2.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h/h)*tmp2/Polynominal::Factorial(i)*dffTab.GetTableValueBackward(i);

#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsGaussForwardInterpolation_ 0
Polynominal Polynominal::SetAsGaussForwardInterpolation(double fn[], int n_set){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAsGaussForwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralForward(0));
	Polynominal tmp(0.0,1.0);
	Pn = Pn + tmp*dffTab.GetTableValueCentralForward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}
#if _DEBUG_Polynominal_SetAsGaussForwardInterpolation_ > 1
		tmp.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + tmp*dffTab.GetTableValueCentralForward(i)/Polynominal::Factorial(i);

#if _DEBUG_Polynominal_SetAsGaussForwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussForwardInterpolation_ 0
Polynominal Polynominal::SetAs1stOrderDerivertiveByGaussForwardInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussForwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralForward(0));
	Pn = Pn.GetDerivertive();

	Polynominal tmp(0.0,1.0);
	Polynominal tmp1 = tmp.GetDerivertive();

	Pn = Pn + (1.0/h)*tmp1*dffTab.GetTableValueCentralForward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}
		Polynominal tmp1 = tmp.GetDerivertive();
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussForwardInterpolation_ > 1
		cout << "tmp : " << tmp.Show('x',-1);
		cout << "tmp1 : " << tmp1.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h)*tmp1/Polynominal::Factorial(i)*dffTab.GetTableValueCentralForward(i);

#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussForwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussForwardInterpolation_ 0
Polynominal Polynominal::SetAs2ndOrderDerivertiveByGaussForwardInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussForwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralForward(0));
	Pn = Pn.GetDerivertive();
	Pn = Pn.GetDerivertive();

	Polynominal tmp(0.0,1.0);
	Polynominal tmp1 = tmp.GetDerivertive();
	Polynominal tmp2 = tmp1.GetDerivertive();

	Pn = Pn + (1.0/h/h)*tmp2*dffTab.GetTableValueCentralForward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}
		Polynominal tmp1 = tmp.GetDerivertive();
		Polynominal tmp2 = tmp1.GetDerivertive();
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussForwardInterpolation_ > 1
		cout << "tmp : " << tmp.Show('x',-1);
		cout << "tmp1 : " << tmp1.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralForward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h/h)*tmp2/Polynominal::Factorial(i)*dffTab.GetTableValueCentralForward(i);

#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussForwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsGaussBackwardInterpolation_ 0
Polynominal Polynominal::SetAsGaussBackwardInterpolation(double fn[], int n_set){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAsGaussBackwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralBackward(0));
	Polynominal tmp(0.0,1.0);
	Pn = Pn + tmp*dffTab.GetTableValueCentralBackward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}
#if _DEBUG_Polynominal_SetAsGaussBackwardInterpolation_ > 1
		tmp.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + tmp*dffTab.GetTableValueCentralBackward(i)/Polynominal::Factorial(i);

#if _DEBUG_Polynominal_SetAsGaussBackwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussBackwardInterpolation_ 0
Polynominal Polynominal::SetAs1stOrderDerivertiveByGaussBackwardInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussBackwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralBackward(0));
	Pn = Pn.GetDerivertive();

	Polynominal tmp(0.0,1.0);
	Polynominal tmp1 = tmp.GetDerivertive();

	Pn = Pn + (1.0/h)*tmp1*dffTab.GetTableValueCentralBackward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}
		tmp1 = tmp.GetDerivertive();
#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussBackwardInterpolation_ > 1
		tmp.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h)*tmp1/Polynominal::Factorial(i)*dffTab.GetTableValueCentralBackward(i);

#if _DEBUG_Polynominal_SetAs1stOrderDerivertiveByGaussBackwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussBackwardInterpolation_ 0
Polynominal Polynominal::SetAs2ndOrderDerivertiveByGaussBackwardInterpolation(double fn[], int n_set, double h){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussBackwardInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(dffTab.GetTableValueCentralBackward(0));
	Pn = Pn.GetDerivertive();
	Pn = Pn.GetDerivertive();

	Polynominal tmp(0.0,1.0);
	Polynominal tmp1 = tmp.GetDerivertive();
	Polynominal tmp2 = tmp1.GetDerivertive();

	Pn = Pn + (1.0/h/h)*tmp2*dffTab.GetTableValueCentralBackward(1);
	for (int i = 2; i < n_set; ++i) {
		if((i % 2)  == 0){
			tmp = tmp*Polynominal((double)(i/2),1.0);
		}else{
			tmp = tmp*Polynominal((double)(-1*i/2),1.0);
		}
		tmp1 = tmp.GetDerivertive();
		tmp2 = tmp1.GetDerivertive();
#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussBackwardInterpolation_ > 1
		tmp.Show('x',-1);
		cout << "Table Value = " << dffTab.GetTableValueCentralBackward(i) << ", Factorial = " << Polynominal::Factorial(i) << endl;
#endif
		Pn = Pn + (1.0/h/h)*tmp2/Polynominal::Factorial(i)*dffTab.GetTableValueCentralBackward(i);

#if _DEBUG_Polynominal_SetAs2ndOrderDerivertiveByGaussBackwardInterpolation_ > 0
		Pn.Show('x',-1);
#endif
	}

	return Pn;
}

#define _DEBUG_Polynominal_SetAsStirlingInterpolation_ 0
Polynominal Polynominal::SetAsStirlingInterpolation(double fn[], int n_set){

	// Get Difference Table at constant difference
	DividedDifferencelTable dffTab(fn, n_set);
#if _DEBUG_Polynominal_SetAsStirlingInterpolation_ > 1
	dffTab.ShowAllTable();
#endif

	Polynominal Pn(0.0);

	if((n_set % 2) == 0){
		cout << "Error : Even number of Points is not implemented.";
	}else{
		Pn = Polynominal(dffTab.GetTableValueCentralForward(0));
		Polynominal tmp_odd(0.0,1.0);
		Pn = Pn + tmp_odd*(dffTab.GetTableValueCentralForward(1)+dffTab.GetTableValueCentralBackward(1))/2.0;
	#if _DEBUG_Polynominal_SetAsStirlingInterpolation_ > 0
		tmp_odd.Show('x',-1);
		Pn.Show('x',-1);
	#endif

		Polynominal tmp_even = tmp_odd.Power(2);
		Pn = Pn + tmp_even*dffTab.GetTableValueCentralForward(2)/Polynominal::Factorial(2);
	#if _DEBUG_Polynominal_SetAsStirlingInterpolation_ > 0
		tmp_even.Show('x',-1);
		Pn.Show('x',-1);
	#endif

		for (int i = 3; i < n_set; ++i) {
			if((i % 2)  == 0){
				tmp_even = tmp_even*Polynominal((double)(i/2),1.0)*Polynominal((double)(-i/2),1.0);
				Pn = Pn + tmp_even*dffTab.GetTableValueCentralForward(i)/Polynominal::Factorial(i);
			#if _DEBUG_Polynominal_SetAsStirlingInterpolation_ > 0
				tmp_even.Show('x',-1);
				Pn.Show('x',-1);
			#endif
			}else{
				tmp_odd = tmp_odd*Polynominal((double)(i/2),1.0)*Polynominal((double)(-i/2),1.0);
				Pn = Pn + tmp_odd*(dffTab.GetTableValueCentralForward(i)+dffTab.GetTableValueCentralBackward(i))/2.0/Polynominal::Factorial(i);
			#if _DEBUG_Polynominal_SetAsStirlingInterpolation_ > 0
				tmp_odd.Show('x',-1);
				Pn.Show('x',-1);
			#endif
			}
		}
	}

	return Pn;
}

// n : Order of Polynominal
Polynominal Polynominal::SetAsLegendre(int n){

	Polynominal Pn;
	Polynominal Pn0(1.0);
	Polynominal Pn1(0.0, 1.0);
	if(n==0){
		Pn = Pn0;
	}else if(n==1){
		Pn = Pn1;
	}else if(n>1){
		for (int i = 2; i <= n; ++i) {
			Polynominal Pn2 = ( (2.0*(i-1) + 1.0)*Polynominal(0.0,1.0)*Pn1 - (i-1)*Pn0 ) / i;
			Pn = Pn2;

			Pn0 = Pn1;
			Pn1 = Pn2;
		}
	}

	return Pn;
}

Polynominal Polynominal::SetAsHermite(int n){

	Polynominal Pn;
	Polynominal Pn0(1.0);
	Polynominal Pn1(0.0, 2.0);
	if(n==0){
		Pn = Pn0;
	}else if(n==1){
		Pn = Pn1;
	}else if(n>1){
		for (int i = 2; i <= n; ++i) {
			Polynominal Pn2 = 2.0*Polynominal(0.0,1.0)*Pn1 - 2.0*(i-1)*Pn0  ;
			Pn = Pn2;

			Pn0 = Pn1;
			Pn1 = Pn2;
		}
	}

	return Pn;
}

Polynominal Polynominal::SetAsLaguerre(int n){

	Polynominal Pn;
	Polynominal Pn0(1.0);
	Polynominal Pn1(1.0, -1.0);
	if(n==0){
		Pn = Pn0;
	}else if(n==1){
		Pn = Pn1;
	}else if(n>1){
		for (int i = 2; i <= n; ++i) {
			Polynominal Pn2 = (2*(i-1) + 1.0 - Polynominal(0.0,1.0))*Pn1 - (i-1)*(i-1)*Pn0  ;
			Pn = Pn2;

			Pn0 = Pn1;
			Pn1 = Pn2;
		}
	}

	return Pn;
}


Polynominal& Polynominal::operator=(const Polynominal& Pn){
	if(this != &Pn){
		this->n = Pn.n;
		delete[] this->an;
		this->an = new double[this->n];
		for (int i = 0; i < this->n; ++i) {
			this->an[i] = Pn.an[i];
		}
	}
	return *this;
}


Polynominal Polynominal::operator +(Polynominal Pn2){


	int nn = Pn2.n > this->n ? Pn2.n : this->n;
	Polynominal ret(nn);

	if(Pn2.n > this->n){
		for (int i = 0; i < this->n; ++i) {
			ret.an[i] = this->an[i] + Pn2.an[i];
		}
		for (int i = this->n; i < Pn2.n; ++i) {
			ret.an[i] = Pn2.an[i];
		}
	}else{
		for (int i = 0; i < Pn2.n; ++i) {
			ret.an[i] = this->an[i] + Pn2.an[i];
		}
		for (int i = Pn2.n; i < this->n; ++i) {
			ret.an[i] = this->an[i];
		}
	}

	return ret.CheckMaxOrder();
}


Polynominal operator +(Polynominal Pn, double a0){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = Pn.an[i];
	}
	ret.an[0] += a0;

	return ret.CheckMaxOrder();
}


Polynominal operator +(double a0, Polynominal Pn){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = Pn.an[i];
	}
	ret.an[0] += a0;

	return ret.CheckMaxOrder();
}


Polynominal Polynominal::operator -(Polynominal Pn2){

	int nn = Pn2.n > this->n ? Pn2.n : this->n;
	Polynominal ret(nn);

	if(Pn2.n > this->n){
		for (int i = 0; i < this->n; ++i) {
			ret.an[i] = this->an[i] - Pn2.an[i];
		}
		for (int i = this->n; i < Pn2.n; ++i) {
			ret.an[i] = -1.0*Pn2.an[i];
		}
	}else{
		for (int i = 0; i < Pn2.n; ++i) {
			ret.an[i] = this->an[i] - Pn2.an[i];
		}
		for (int i = Pn2.n; i < this->n; ++i) {
			ret.an[i] = this->an[i];
		}
	}

	return ret.CheckMaxOrder();
}


Polynominal operator -(Polynominal Pn, double a0){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = Pn.an[i];
	}
	ret.an[0] -= a0;

	return ret.CheckMaxOrder();
}


Polynominal operator -(double a0, Polynominal Pn){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = -1.0*Pn.an[i];
	}
	ret.an[0] += a0;

	return ret.CheckMaxOrder();
}


Polynominal Polynominal::operator *(Polynominal Pn2){

	int nn;
	if((Pn2.n == 1) || (this->n == 1)){
		 nn = Pn2.n > this->n ? Pn2.n : this->n;
	}else{
		nn = (Pn2.n-1)+(this->n-1)+1;
	}

	Polynominal ret(nn);

	for (int i = 0; i < this->n; ++i) {
		for (int j = 0; j < Pn2.n; ++j) {
			ret.an[i+j] += this->an[i]*Pn2.an[j];
		}
	}

	return ret;
}


Polynominal operator *(Polynominal Pn, double a0){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = Pn.an[i]*a0;
	}

	return ret;
}


Polynominal operator *(double a0, Polynominal Pn){

	int nn = Pn.n;
	Polynominal ret(nn);

	for (int i = 0; i < nn; ++i) {
		ret.an[i] = a0*Pn.an[i];
	}

	return ret;
}

// Get PnNum/PnDen(n@PnNum > n@PnDen), PnNum = PnDen*PnRet + PnRest
//Polynominal SyntheticDivision(Polynominal PnNum, Polynominal PnDen, Polynominal *PnRest){
//}


Polynominal operator/(Polynominal Pn, double a){
	Polynominal ret(Pn.n);

	for (int i = 0; i < Pn.n; ++i) {
		ret.an[i] = Pn.an[i]/a;
	}

	return ret;
}


Polynominal Polynominal::Power(unsigned int pow){

	Polynominal ret(1.0);

	if(pow > 0){
		for (unsigned int i = 0; i < pow; ++i) {
			ret = ret * (*this);
		}
	}

	return ret.CheckMaxOrder();
}


Polynominal Polynominal::Substitute(Polynominal Pn){

	Polynominal ret(0.0);

	for (int i = 0; i < this->n; ++i) {
		ret = ret + this->an[i]*Pn.Power((unsigned int)i);
	}

	return ret.CheckMaxOrder();
}


Polynominal Polynominal::GetDerivertive(){

	Polynominal ret;

	if(this->n > 1){
		int nn = this->n-1;
		 ret =Polynominal(nn);

		for (int i = 0; i < nn; ++i) {
			ret.an[i] = this->an[i+1]*(i+1);
		}
	}else{
		 ret =Polynominal(0.0);
	}

	return ret;
}


Polynominal Polynominal::GetIntegral(){

	Polynominal ret;

	int nn = this->n+1;
	ret =Polynominal(nn);

	ret.an[0] = 0.0;
	for (int i = 1; i < nn; ++i) {
		ret.an[i] = this->an[i-1]/(double)i;
	}

	return ret;
}



// Get value of f(x) for given x
double Polynominal::GetValue(double x){
#define _Polynominal_GetValue_ 1
#if _Polynominal_GetValue_ == 0
	double f=0.0;
	for(int i=0; i<this->n; i++){
		f += this->an[i]*pow(x,i);
	}
#elif _Polynominal_GetValue_ == 1
	double f = this->an[this->n-1];
	for(int i=this->n-2; i>=0; i--){
		f = f*x + this->an[i];
	}
#endif


	return f;
}

int Polynominal::Factorial(int n){

	int ret = 1;

	for (int i = 0; i < n; ++i) {
		ret *= n-i;
	}

	return ret;
}

void Test_Class_Polynominal(){

#define _DEBUG_Test_Class_Polynominal_ 14
#if _DEBUG_Test_Class_Polynominal_ == 1
	Polynominal Pn1 =Polynominal(1.0); Pn1.Show('x',-1);
	Polynominal Pn2 = Polynominal(2.0,1.0); Pn2.Show('x',-1);
	Polynominal Pn3 = Pn1+Pn2; Pn3.Show('x',-1);
	Polynominal Pn4 = Pn2-Pn1; Pn4.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 2
	Polynominal Pn1 =Polynominal(1.0,1.0); Pn1.Show('x',-1);
	Polynominal Pn2 = Polynominal(1.0,1.0); Pn2.Show('x',-1);
	Polynominal Pn3 = Pn1*Pn2; Pn3.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 3
	Polynominal Pn1 =Polynominal(2.0,4.0); Pn1.Show('x',-1);
	Polynominal Pn2 = Pn1 / 2.0; Pn2.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 4
// Test Lagrange Interpolate
	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	Polynominal Pn = Polynominal::SetAsLagrangeInterpolation(xn,fn,n); Pn.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 5
// Test Get Derivertive function
	int n=4;
	double an[] = {0.0, 1.0, 2.0, 3.0};
	Polynominal Pn1 =Polynominal(n, an); Pn1.Show('x',-1);
	Polynominal Pn2 = Pn1.GetDerivertive(); Pn2.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 6
// Test Hermite Interpolate
	int n=2;
	double xn[] = {0.0, 1.0};
	double fn[] = {0.0, 2.0};
	double fpn[] = {1.0, 3.0};
	Polynominal Pn = Polynominal::SetAsHermiteInterpolation(xn, fn, fpn, n); Pn.Show('x',-1);
//	Pn.CheckMaxOrder().Show('x',-1);
//	Polynominal::SetAsHermiteInterpolation(xn, fn, fpn, n).Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 7
// Test Newton Interpolate
	int n=4;
	double xn[] = {0.0, 1.0, 2.0, 3.0};
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	Polynominal Pn = Polynominal::SetAsNewtonInterpolation(xn,fn,n); Pn.Show('x',-1);
//	Polynominal::SetAsNewtonInterpolation(xn,fn,n).Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 8
// Test Newton Forward Difference Interpolation and get derivertive
	int n=4;
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	Polynominal Pn = Polynominal::SetAsNewtonForwardDifferenceInterpolation(fn,n); Pn.Show('x',-1);
//	Polynominal::SetAsNewtonForwardDifferenceInterpolation(fn,n).Show('x',-1);

	cout << "Derivertive" << endl;
	Polynominal Pn1 = Pn.GetDerivertive(); Pn1.Show('x',-1);
	cout << "Pn'(2.5) = " << Pn1.GetValue(2.5) << endl;
	Polynominal Pn12 = Polynominal::SetAs1stOrderDerivertiveByNewtonForwardDifferenceInterpolation(fn,n,1.0); Pn12.Show('x',-1);
	cout << "Pn'(2.5) = " << Pn12.GetValue(2.5) << endl;

	Polynominal Pn2 = Pn1.GetDerivertive(); Pn2.Show('x',-1);
	cout << "Pn''(2.5) = " << Pn2.GetValue(2.5) << endl;
	Polynominal Pn22 = Polynominal::SetAs2ndOrderDerivertiveByNewtonForwardDifferenceInterpolation(fn,n,1.0); Pn22.Show('x',-1);
	cout << "Pn''(2.5) = " << Pn22.GetValue(2.5) << endl;

	cout << "Integral" << endl;
	Polynominal IntegPn = Pn.GetIntegral(); IntegPn.Show('x',-1);

#elif _DEBUG_Test_Class_Polynominal_ == 9
	// Test Newton Backward Difference Interpolation
	int n=4;
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	Polynominal Pn = Polynominal::SetAsNewtonBackwardDifferenceInterpolation(fn,n); Pn.Show('x',-1);
	Polynominal Pn01(-3.0,1.0); Pn01.Show('x',-1);
	Polynominal Pn02 = Pn.Substitute(Pn01); Pn02.Show('x',-1);

	Polynominal Pn1 = Pn.GetDerivertive(); Pn1.Show('x',-1);
	cout << "Pn'(-0.5) = " << Pn1.GetValue(-0.5) << endl;
	Polynominal Pn12 = Polynominal::SetAs1stOrderDerivertiveByNewtonBackwardDifferenceInterpolation(fn,n,1.0); Pn12.Show('x',-1);
	cout << "Pn'(-0.5) = " << Pn12.GetValue(-0.5) << endl;

	Polynominal Pn2 = Pn1.GetDerivertive(); Pn2.Show('x',-1);
	cout << "Pn''(-0.5) = " << Pn2.GetValue(-0.5) << endl;
	Polynominal Pn22 = Polynominal::SetAs2ndOrderDerivertiveByNewtonBackwardDifferenceInterpolation(fn,n,1.0); Pn22.Show('x',-1);
	cout << "Pn'(-0.5) = " << Pn22.GetValue(-0.5) << endl;
#elif _DEBUG_Test_Class_Polynominal_ == 10
// Test Power function
	Polynominal Pn1(1.0,1.0); Pn1.Show('x',-1);
	Polynominal Pn2 = Pn1.Power(0); Pn2.Show('x',-1);
	Pn2 = Pn1.Power(1); Pn2.Show('x',-1);
	Pn2 = Pn1.Power(2); Pn2.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 11
// Test Substitute function
	int n=4;
	double fn[] = {1.0, 4.0, 15.0, 40.0};
	Polynominal Pn1 = Polynominal::SetAsNewtonForwardDifferenceInterpolation(fn,n); Pn1.Show('x',-1);
	Polynominal Pn2 = Pn1.Substitute(Polynominal(1.0)); Pn2.Show('x',-1);
	Polynominal Pn31(0.0,1.0); Pn31.Show('x',-1);
	Polynominal Pn32 = Pn1.Substitute(Pn31); Pn32.Show('x',-1);
	int nn=3;
	double an[] = {0.0, 0.0, 1.0};
	Polynominal Pn41(nn,an); Pn41.Show('x',-1);
	Polynominal Pn42 = Pn1.Substitute(Pn41); Pn42.Show('x',-1);
#elif _DEBUG_Test_Class_Polynominal_ == 12
// Test Gauss Forward/Backward Interpolation
	int n=4;
//	double fn[] = {1.0, 4.0, 15.0, 40.0};
	double fn[] = {sin(0.0), sin(0.1), sin(0.2), sin(0.3)};

	cout << "Gauss Forward Interpolation" << endl;
	Polynominal Pn1 = Polynominal::SetAsGaussForwardInterpolation(fn,n); Pn1.Show('x',-1);
	cout << Pn1.GetValue(0.5) << endl;

	cout << "Gauss Backward Interpolation" << endl;
	Polynominal Pn2 = Polynominal::SetAsGaussBackwardInterpolation(fn,n); Pn2.Show('x',-1);
	cout << Pn2.GetValue(-0.5) << endl;
#elif _DEBUG_Test_Class_Polynominal_ == 13
// Test Stirling Interpolation
	int n=5;
	double fn[] = {sin(-0.1), sin(0.0), sin(0.1), sin(0.2), sin(0.3)};

	cout << "Stirling Interpolation" << endl;
	Polynominal Pn = Polynominal::SetAsStirlingInterpolation(fn,n); Pn.Show('x',-1);
	cout << Pn.GetValue(0.5) << endl;
#elif _DEBUG_Test_Class_Polynominal_ == 14
	// Test Get Derivertive by Gauss Forward/Backward Interpolation
		int n=4;
		double fn[] = {1.0, 4.0, 15.0, 40.0};

		cout << "Gauss Forward Interpolation" << endl;
		Polynominal Pn0f = Polynominal::SetAsGaussForwardInterpolation(fn,n); Pn0f.Show('x',-1);
		cout << "Pn(1.5)=" << Pn0f.GetValue(1.5) << endl;
		Polynominal Pn1f1 = Pn0f.GetDerivertive(); Pn1f1.Show('x',-1);
		cout << "Pn'(1.5)=" << Pn1f1.GetValue(1.5) << endl;
		Polynominal Pn1f2 = Polynominal::SetAs1stOrderDerivertiveByGaussForwardInterpolation(fn,n,1.0); Pn1f2.Show('x',-1);
		cout << "Pn'(1.5)=" << Pn1f2.GetValue(1.5) << endl;
		Polynominal Pn2f1 = Pn1f1.GetDerivertive(); Pn2f1.Show('x',-1);
		cout << "Pn''(1.5)=" << Pn2f1.GetValue(1.5) << endl;
		Polynominal Pn2f2 = Polynominal::SetAs2ndOrderDerivertiveByGaussForwardInterpolation(fn,n,1.0); Pn2f2.Show('x',-1);
		cout << "Pn''(1.5)=" << Pn2f2.GetValue(1.5) << endl;

		cout << "Gauss Backward Interpolation" << endl;
		Polynominal Pn0b = Polynominal::SetAsGaussBackwardInterpolation(fn,n); Pn0b.Show('x',-1);
		cout << "Pn(0.5)=" << Pn0b.GetValue(0.5) << endl;
		Polynominal Pn1b1 = Pn0b.GetDerivertive(); Pn1b1.Show('x',-1);
		cout << "Pn'(0.5)=" << Pn1b1.GetValue(0.5) << endl;
		Polynominal Pn1b2 = Polynominal::SetAs1stOrderDerivertiveByGaussBackwardInterpolation(fn,n,1.0); Pn1b2.Show('x',-1);
		cout << "Pn'(0.5)=" << Pn1b2.GetValue(0.5) << endl;
		Polynominal Pn2b1 = Pn1b1.GetDerivertive(); Pn2b1.Show('x',-1);
		cout << "Pn''(0.5)=" << Pn2b1.GetValue(0.5) << endl;
		Polynominal Pn2b2 = Polynominal::SetAs2ndOrderDerivertiveByGaussBackwardInterpolation(fn,n,1.0); Pn2b2.Show('x',-1);
		cout << "Pn''(0.5)=" << Pn2b2.GetValue(0.5) << endl;
#elif _DEBUG_Test_Class_Polynominal_ == 15
// Test SetAsLegendre
		int n = 5;
		Polynominal Pn[n];
		for (int i = 0; i <=n; ++i) {
			cout << i << endl;
			Pn[i] = Polynominal::SetAsLegendre(i); Pn[i].Show('x',-1);
		}
#elif _DEBUG_Test_Class_Polynominal_ == 16
// Test SetAsHetmite
		int n = 5;
		Polynominal Pn[n];
		for (int i = 0; i <=n; ++i) {
			cout << i << endl;
			Pn[i] = Polynominal::SetAsHermite(i); Pn[i].Show('x',-1);
		}
#elif _DEBUG_Test_Class_Polynominal_ == 17
// Test SetAsLaguerre
		int n = 5;
		Polynominal Pn[n];
		for (int i = 0; i <=n; ++i) {
			cout << i << endl;
			Pn[i] = Polynominal::SetAsLaguerre(i); Pn[i].Show('x',-1);
		}
#endif
	return;
}


class GaussIntegralFormula{
private:
	int n;
	double *xk;
	double *wk;

public:
	GaussIntegralFormula();
	GaussIntegralFormula(const GaussIntegralFormula&); // Copy Constructor for Initialize
	GaussIntegralFormula(int n_set, double xk_set[], double wk_set[]);
	GaussIntegralFormula(int NumPoly);
	~GaussIntegralFormula(){
		if(this->xk != NULL){ delete[] this->xk; }
		if(this->wk != NULL){ delete[] this->wk; }
		this->n=0;
	};
	void Show(char);
	int GetDimenstion(){ return this->n; }
	double GetLegendreXk(int nplus1, int idx);
};

GaussIntegralFormula::GaussIntegralFormula(){
	this->n = 0;
	this->xk = NULL;
	this->wk = NULL;
}

GaussIntegralFormula::GaussIntegralFormula(const GaussIntegralFormula &gIntg){
	this->n = gIntg.n;
	this->xk = new double[this->n];
	this->wk = new double[this->n];
	for (int i = 0; i < this->n; ++i) {
		this->xk[i] = gIntg.xk[i];
		this->wk[i] = gIntg.wk[i];
	}
}

GaussIntegralFormula::GaussIntegralFormula(int n_set, double xk_set[], double wk_set[]){
	this->n = n_set;
	this->xk = new double[this->n];
	this->wk = new double[this->n];

	for (int i = 0; i < this->n; ++i) {
		this->xk[i] = xk_set[i];
		this->wk[i] = wk_set[i];
	}

}

#define _File_GaussIntegralFormula_Legendre_ "C:\\Users\\Kani\\Work\\Introduction and Applications to Numerical Calculation by C++\\src\\GaussIntegralFormula_Legendre.csv"
#define _File_GaussIntegralFormula_Hermite_ "C:\\Users\\Kani\\Work\\Introduction and Applications to Numerical Calculation by C++\\src\\GaussIntegralFormula_Hermite.csv"
#define _File_GaussIntegralFormula_Laguerre_ "C:\\Users\\Kani\\Work\\Introduction and Applications to Numerical Calculation by C++\\src\\GaussIntegralFormula_Laguerre.csv"

// NumPoly: 1=Legendre, 2=Hermite, 3=Laguerre
GaussIntegralFormula::GaussIntegralFormula(int NumPoly){
//	char *FileName;
	string FileName;

	switch(NumPoly){
	case 1:
		FileName = _File_GaussIntegralFormula_Legendre_;
		break;
	case 2:
		FileName = _File_GaussIntegralFormula_Hermite_;
		break;
	case 3:
		FileName = _File_GaussIntegralFormula_Laguerre_;
		break;
	}
	ifstream fin(FileName);
	if(!fin){
		cout << "Couldn't Open file \"" << FileName << "\"" << endl;
		return ;
	}

    string str;
    while(getline(fin,str)){
        string token;
        istringstream stream(str);

        //1行のうち、文字列とコンマを分割する
        while(getline(stream,token,',')){
            //すべて文字列として読み込まれるため
            //数値は変換が必要
            double temp=std::stod(token);
            cout<<temp<<",";
        }
        cout<<endl;
    }
	fin.close();
}

void GaussIntegralFormula::Show(char delimiter = ','){

	cout << "xk : ";
	for(int i=0; i<this->n; i++){
		cout << this->xk[i];
		if(i<this->n-1){ cout << delimiter; }
	}
	cout << endl;

	cout << "wk : ";
	for(int i=0; i<this->n; i++){
		cout << this->wk[i];
		if(i<this->n-1){ cout << delimiter; }
	}
	cout << endl;

}

// nplus1:n+1=2,3,4,5,6, i=1,2,...,nplus1
double GaussIntegralFormula::GetLegendreXk(int nplus1, int idx){

	int start=0;
	for (int i = 2; i <= nplus1; ++i) {
		start+=i;
	}
	cout << start << endl;

	return this->xk[start-1+idx-1];
}


#define _DEBUG_Test_Class_GaussIntegralFormula_ 1
void Test_Class_GaussIntegralFormula(){

	GaussIntegralFormula gIntg=GaussIntegralFormula(1);

//	cout << "n=" << gIntg.GetDimenstion() << endl;
//	gIntg.Show();

//	cout << gIntg.GetLegendreXk(2,1) << endl;


#if _DEBUG_Test_GaussIntegralFormula_ == 1
#endif

	return;
}

