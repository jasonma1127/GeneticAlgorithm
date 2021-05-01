#include <iostream>
#include <random>
#include <ctime>
#include <algorithm>
#include <fstream>

using namespace std;

//Problem requires
int maxValue_x = 10, minValue_x = 8;
int maxValue_y = 13, minValue_y = 10;
int constrain = 22;
int precision = 3;

//Algorithm variables
int generation = 3000;
int population_size = 350;
double crossoverRate = 0.95;
double mutationRate = 0.2;

//輸出文字檔 
ofstream newFile;	

double function(double x, double y){
	return -x * sin(4 * x)-1.1*y*sin(2*y)+1;
}

int xbits, ybits;	//x與y之bits的長度 
int chrom_length;

int **arr = new int*[population_size];// 宣告二維陣
default_random_engine generator(time(NULL));//亂數產生器
uniform_real_distribution<float> unif(0.0, 1.0);// 亂數的機率分布
vector<double> fitnessVector;

double localOptAns = -100;
double localOpt_x, localOpt_y;
double globalOptAns = -100;
double globalOpt_x, globalOpt_y;

vector<int> parentVector_1;
vector<int> parentVector_2;

int binaryToDecimal(vector<int> _tmpVector){
	int num;
	reverse(_tmpVector.begin(), _tmpVector.end()); 
	
	for(int i = 0; i < _tmpVector.size(); i ++){
		num = num + _tmpVector[i] * pow(2, i);
	}
	return num;
}

int checkConstrain(double _realNum_x, double _realNum_y, int numOfChrom){
	double funcAns;
	if(_realNum_x + _realNum_y > constrain){
		//cout<<"constrain alert!!"<<endl;
		for (int i = 0; i < chrom_length; i++)
		{
			arr[numOfChrom][i] = round(unif(generator));// 產生亂數
		}
		return numOfChrom-1;		
	}else{
		funcAns = function(_realNum_x, _realNum_y);
		fitnessVector.push_back(funcAns);
		if(localOptAns < funcAns){
			localOptAns = funcAns;
			localOpt_x = _realNum_x;
			localOpt_y = _realNum_y;
		}
		return numOfChrom;
	}
}

void setUp(){
	
	//計算X長度及需要的bits
	for (int i = 0; ; i++)
	{
		if (pow(2, i) > ((maxValue_x - minValue_x)*pow(10, precision))) {
			xbits = i;
			break;
		}
	}

	//計算Y長度及需要的bits
	for (int i = 0; ; i++)
	{
		if (pow(2, i) > ((maxValue_y - minValue_y)*pow(10, precision))) {
			ybits = i;
			break;
		}
	}
	
	//======================//
	xbits = 12; ybits = 12;
	
	chrom_length = xbits + ybits;
}

void initialization(){

	cout<<"init: "<<xbits<<" + "<<ybits<<" = "<<chrom_length<<endl;
	
	for (int i = 0; i < population_size; i++) {
		arr[i] = new int[chrom_length];
	}

	for (int i = 0; i < population_size; i++) {
		for (int j = 0; j < chrom_length; j++) {
			
			arr[i][j] = round(unif(generator));// 產生亂數
		}
	}
}

void evalution(){
		
	int decimal_x, decimal_y;
	double realNum_x, realNum_y;
	vector<int> tmpVector;
	
	for(int i = 0; i < population_size; i ++){
		
		tmpVector.clear();
		for(int j = 0; j < xbits; j ++){
			tmpVector.push_back(arr[i][j]);
		}
		
		decimal_x = binaryToDecimal(tmpVector);
		realNum_x = (minValue_x - 0.5) +  decimal_x * ((maxValue_x - minValue_x) / ((pow(2, xbits) - 1))); //輸出realNum_x值
		
		tmpVector.clear();
		for(int j = ybits; j < chrom_length; j ++){
			tmpVector.push_back(arr[i][j]);
		}
		
		decimal_y = binaryToDecimal(tmpVector);
		realNum_y = (minValue_y - 0.5) +  decimal_y * ((maxValue_y - minValue_y) / ((pow(2, ybits) - 1))); //輸出realNum_x值
		
		i = checkConstrain(realNum_x, realNum_y, i);		
	}

	if(globalOptAns < localOptAns){
		globalOptAns = localOptAns;
		globalOpt_x = localOpt_x;
		globalOpt_y = localOpt_y;
	}
	cout<<"LocalOpt_x: "<<localOpt_x<<" LocalOpt_y: "<<localOpt_y<<" LocalOptAns: "<<localOptAns<<endl;
	newFile << localOptAns << " " << globalOptAns << endl;
}

void selection(){
	
	int parent_1, parent_2;
		
	//parent_1
	//亂數選出select_1與select_2來競爭
	double rand = unif(generator);
	int select_1 = floor(rand*(population_size));
	
	rand = unif(generator);
	int select_2 = floor(rand*(population_size));
	
	//當select_1與select_2相同則重選到不同 
	while(select_1 == select_2){
		rand = unif(generator);
		select_2 = floor(rand*(population_size));
	}
	
	//兩者競爭，數值大者設為parent_1 
	if(fitnessVector[select_1] > fitnessVector[select_2]){
		parent_1 = select_1;
	}else{
		parent_1 = select_2;
	}
	
	//在Initialization中的群體找尋parent_1的染色體 
	parentVector_1.clear();
	for (int j = 0; j < chrom_length; j++)
	{
		parentVector_1.push_back(arr[parent_1][j]);
	}
	
	//把parent_1的染色體從群體中刪除 
	for (int i = parent_1; i < population_size-1; i++)
	{
		for (int j = 0; j < chrom_length; j++)
		{
			arr[i][j] = arr[i + 1][j];
		}
	}
	
	//parent_2	
	//亂數選出select_1與select_2來競爭 
	rand = unif(generator);
	select_1 = floor(rand*(population_size - 1));
	
	rand = unif(generator);
	select_2 = floor(rand*(population_size - 1));
	
	//當select_1與select_2相同則重選到不同 
	while(select_1 == select_2){	
		rand = unif(generator);
		select_2 = floor(rand*(population_size - 1));
	}
	
	//兩者競爭，數值大者設為parent_2
	if(fitnessVector[select_1] > fitnessVector[select_2]){
		parent_2 = select_1;
	}else{
		parent_2 = select_2;
	}
	
	//在Initialization中的群體找尋parent_2的染色體 
	parentVector_2.clear();
	for (int j = 0; j < chrom_length; j++)
	{
		parentVector_2.push_back(arr[parent_2][j]);
	}
	
	//把parent_2的染色體從群體中刪除
	for (int i = parent_2; i < population_size-1; i++)
	{
		for (int j = 0; j < chrom_length; j++)
		{
			arr[i][j] = arr[i + 1][j];
		}
	}
}

void crossover(){
	
	double cutlocation_randnum;

	//亂數產生x的cutPoint與y的cutPoint
	cutlocation_randnum = unif(generator);
	int cut_location_x = ceil(cutlocation_randnum*(xbits - 1));
	
	cutlocation_randnum = unif(generator);
	int cut_location_y = ceil(cutlocation_randnum*(ybits - 1)) + xbits;
	
	//將要交配的x子染色體序列與y子染色體序列取出 
	vector<int> tmp_subvector_x;
	tmp_subvector_x = vector<int>(parentVector_1.begin() + cut_location_x, parentVector_1.begin() + xbits);	
	
	vector<int> tmp_subvector_y;
	tmp_subvector_y = vector<int>(parentVector_1.begin() + cut_location_y, parentVector_1.end());

	//x與y個別進行交配
	parentVector_1.erase(parentVector_1.begin() + cut_location_x, parentVector_1.begin() + xbits);
	parentVector_1.insert(parentVector_1.begin() + cut_location_x, parentVector_2.begin() + cut_location_x, parentVector_2.begin() + xbits);
				
	parentVector_1.erase(parentVector_1.begin() + cut_location_y, parentVector_1.end());
	parentVector_1.insert(parentVector_1.end(), parentVector_2.begin() + cut_location_y, parentVector_2.end());
	
	parentVector_2.erase(parentVector_2.begin() + cut_location_x, parentVector_2.begin() + xbits);
	parentVector_2.insert(parentVector_2.begin() + cut_location_x, tmp_subvector_x.begin(), tmp_subvector_x.end());
	
	parentVector_2.erase(parentVector_2.begin() + cut_location_y, parentVector_2.end());
	parentVector_2.insert(parentVector_2.end(), tmp_subvector_y.begin(), tmp_subvector_y.end());
}

void mutation(vector<int> *_parentVector){
			
	double mut_location_randnum_1 = unif(generator);
	
	//亂數產生變異的染色體位置 
	int mut_location_1 = ceil(mut_location_randnum_1*chrom_length);

	//若原本為0則變1，原本為1則變0 
	if (_parentVector->at(mut_location_1 - 1) == 0)
	{
		replace(_parentVector->begin() + mut_location_1 - 1, _parentVector->end() - ((chrom_length - 1) - (mut_location_1 - 1)), 0, 1);
	}
	else
	{
		replace(_parentVector->begin() + mut_location_1 - 1, _parentVector->end() - ((chrom_length - 1) - (mut_location_1 - 1)), 1, 0);
	}
}

void replaceParent(){
	//處理完後的染色體存回群體
	for (int i = 0; i < chrom_length; i++)
	{
		arr[population_size - 2][i] = parentVector_1[i];
	}
	for (int i = 0; i < chrom_length; i++)
	{
		arr[population_size - 1][i] = parentVector_2[i];
	}
}

int main(){
	
	setUp();	
	initialization();
	newFile.open("final.txt");

	
	for(int gen = 0; gen < generation; gen ++){
		
		cout<<"==========第 "<<gen+1<<" 次==========="<<endl; 
		
		fitnessVector.clear();
		
		evalution();
		selection();
		
		double doCrossoverRate = unif(generator);
		if(doCrossoverRate <= crossoverRate){
			
			crossover();
			
			double doMutationRate_1 = unif(generator);
			double doMutationRate_2 = unif(generator);
			
			if(doMutationRate_1 <= mutationRate){
				mutation(&parentVector_1);
			}
			if(doMutationRate_2 <= mutationRate){
				mutation(&parentVector_2);	
			}
		}
		replaceParent();
	}
	
	cout<<"=================================================================="<<endl;
	cout<<"GlobalOpt_x: "<<globalOpt_x<<" GlobalOpt_y: "<<globalOpt_y<<" GlobalOptAns: "<<globalOptAns<<endl;
}
