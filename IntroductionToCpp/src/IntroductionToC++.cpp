//============================================================================

// Name        : IntroductionToC++.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstring>
#include <iostream>   // std::cout
#include <string>     // std::string, std::stod

using namespace std;

void Test_NumberSystem(){

	cout << "10 in Binary Number System is " << 0b10 << "\n";
	cout << "10 in Octal Number System is " << 010 << "\n";
	cout << "10 in Decimal Number System is " << 10 << "\n";
	cout << "10 in HexaDecimal Number System is " << 0x10 << "\n";

	return;
}

void Test_SizeOf() {

	int a=1;
	double b=2;

	cout << "size of short integer = " << sizeof(short int) << "byte\n";
	cout << "size of integer = " << sizeof(int) << "byte\n";
	cout << "size of long integer = " << sizeof(long int) << "byte\n";
	cout << "size of float = " << sizeof(float) << "byte\n";
	cout << "size of double = " << sizeof(double) << "byte\n";
	cout << "size of long double = " << sizeof(long double) << "byte\n";

	cout << "size of variable a = " << sizeof(a) << "byte\n";
	cout << "size of variable b = " << sizeof(b) << "byte\n";
	cout << "size of equation a+b = " << sizeof(a+b) << "byte\n";

	return;
}

void Test_ShiftOperator(){

	int a =	2;
	int b,c;

	cout << "a = " << a << "\n";

	b= a<<2;
	cout << "b = " << b << "\n";

	c= a>>1;
	cout << "c = " << c << "\n";

	return;
}

void Test_Cast(){

	int num1=5,num2=4;
	double div;

	div=num1/num2;
	cout << num1 << "/" << num2 << "=" << div <<"\n";

	div=(double)(num1/num2);
	cout << num1 << "/" << num2 << "=" << div <<"\n";

	div=(double)num1/(double)num2;
	cout << num1 << "/" << num2 << "=" << div <<"\n";

	return;
}

void Test_InnerBlockVariable(){

	int i=23;

	for(int i=0; i<50;i++){
		cout << "i = " << i << "\n";
	}

	cout << "i = " << i << "\n";

	return;
}

inline int max(int x,int y){
	if(x>y)
		return x;
	else
		return y;
}

void Test_InlineFunction(){

	int a=101, b=201;

	cout << "a=" << a << ", b=" << b << ", max(a,b)=" << max(a,b) << "\n";

	return;
}

int min(int x,int y){
	if(x<y)
		return x;
	else
		return y;
}

double min(double x,double y){
	if(x<y)
		return x;
	else
		return y;
}

void Test_FunctionOverLoad(){
	int a,b;
	double da,db;

	cout << "Input 2 integer number.";
	cin >> a >> b;
	cout << "a=" << a << ", b=" << b << "\n";

	cout << "Input 2 double number.";
	cin >> da >> db;
	cout << "da=" << da << ", db=" << db << "\n";

	int ans1=min(a,b);
	double ans2=min(da,db);

	cout << "Minimum of integer numbers = " << ans1 << "\n";
	cout << "Minimum of double numbers = " << ans2 << "\n";

	return;

}

template <class T>
T mint(T x, T y){
	if(x<y)
		return x;
	else
		return y;

}

void Test_Template(){
	int a,b;
	double da,db;

	cout << "Input 2 integer number.";
	cin >> a >> b;
	cout << "a=" << a << ", b=" << b << "\n";

	cout << "Input 2 double number.";
	cin >> da >> db;
	cout << "da=" << da << ", db=" << db << "\n";

	//int ans1=min(a,b);
	// double ans2=min(da,db);
	int ans1=mint(a,b);
	double ans2=mint(da,db);

	cout << "Minimum of integer numbers = " << ans1 << "\n";
	cout << "Minimum of double numbers = " << ans2 << "\n";

	return;

}

void Test_Pointer(){
	int a;
	int *pa=&a;

	a=123;

	cout << "Value of variable a is " << a << "\n";
	cout << "Address of variable a is " << &a << "\n";
	cout << "Value of pointer pa is " << pa << "\n";
	cout << "Value of pointer *pa is " << *pa << "\n";

	return;
}

void swap(int *px,int *py){
	int tmp;

	tmp = *px;
	*px = *py;
	*py = tmp;

	return;
}

void swap(int &x,int &y){
	int tmp;

	tmp = x;
	x = y;
	y = tmp;

	return;
}

void Test_Pointer_Swap(){
	int num1 = 5;
	int num2 = 20;

	cout << "Value of variable num1 is " << num1 << "\n";
	cout << "Value of variable num2 is " << num2 << "\n";
	cout << "Swap value of num1 and num2\n";

	//swap(&num1,&num2);
	swap(num1,num2);

	cout << "Value of variable num1 is " << num1 << "\n";
	cout << "Value of variable num2 is " << num2 << "\n";

}

double avg(int t[],int num){
	double sum = 0;

	for(int i=0; i<num; i++){
		sum += t[i];
	}
	return sum/num;
}

void Test_Array_Average(){
	int num = 5;
	int test[num];

	cout << "Input Test Scores of "<< num << " members.\n";
	for(int i=0; i<num; i++){
		cin >> test[i];
	}

	double ans = avg(test,num);
	cout << "Average Test Score is " << ans << "\n";

	return;
}

void Test_String(){
	char str[]="Hello Introduction to C++";
	char *str1;

	cout << str << endl;

	str1="Hello Saltwater Fishing!";

	cout << str1 << endl;

	str1="Aging, Mebaring are my Preasure!!";

	cout << str1 << endl;

	for(int i=0; str1[i]!='\0'; i++){
		cout << str1[i] << '*';
	}
	cout << '\n';

	cout << "Length of string1 is " << strlen(str1) << ".\n";
	cout << "Compare 2 string : " << strcmp(str1,str1) << endl;

	char str2[20], str3[10], str4[10];
	strcpy(str3,"Hello");
	strcpy(str4,"Goodbye");
	cout << str3 << str4 << endl;
	strcpy(str2,str3);
	cout << str2 << endl;
	strcat(str2,str4);
	cout << str2 << endl;

}

int a=0;
void Test_Variable_Scope(){

	int a=0;
	::a=10;

	cout << "local a=" << a << endl;
	cout << "global a=" << ::a << endl;

	a++;
	::a--;

	cout << "local a=" << a << endl;
	cout << "global a=" << ::a << endl;
}

void Func_Variable_LifeTime(){
	int b=0;
	static int c=0;

	cout << "Variable a is " << a << ", Variable b is " << b << ", Variable c is " << c << endl;

	a++;
	b++;
	c++;
}

void Test_Variable_LifeTime(){

	for(int i=0; i<5; i++){
		cout << "i is " << i << endl;
		Func_Variable_LifeTime();
	}
	return;
}



int Test_stod ()
{
  std::string orbits ("365.24 29.53");
  std::string::size_type sz;     // alias of size_t

  double earth = std::stod (orbits,&sz);
  double moon
  = std::stod (orbits.substr(sz));
  std::cout << "The moon completes " << (earth/moon) << " orbits per Earth year.\n";
  return 0;
}

int main(){
	/*
	 * Easy Introduction to C++
	*/

	// Test_SizeOf();
	// Test_ShiftOperator();
	// Test_Cast();
	 Test_NumberSystem();
	// Test_InnerBlockVariable();

	// Test_InlineFunction();
	// Test_FunctionOverLoad();
	// Test_Template();
	// Test_Pointer();
	// Test_Pointer_Swap();
//	 Test_Array_Average();
	//	Test_String();
	//	Test_Scope();
//	 Test_Variable_LifeTime();

//	Test_stod();
}


