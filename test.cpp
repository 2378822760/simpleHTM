#include "encoder.hpp"
#include "SpatialPooler.hpp"
#include "TemporalMemory.hpp"
#include <iostream>
#include <cstdio>



using namespace std;

/**
 * 返回arr的上下界，上界保存在max，下界保存在min
 */
void getArrayRange(const double arr[], int arraySize, double& min, double& max) {
	// 设计函数同时获得最小最大值
	min = max = arr[0];

	for (int i = 0; i < arraySize; ++i) {
		double tmp = arr[i];
		min = min > tmp ? tmp : min;
		max = max < tmp ? tmp : max;
	}

	max += (max - min) / (arraySize - 1);
}

void printCodedDate(UInt inputVector[], UInt vectorSize) {
	for (int i = 0; i < vectorSize; ++i) {
		cout << inputVector[i];
	}
	cout << endl;
}

int main() {
	/* Input */
	double date[] = { 10,10.5,11,11.5,12,12.5,13,14,15,16,17,18,19 };
	int arraySize = sizeof(date) / sizeof(double);
	double minElement, maxElement;
	getArrayRange(date, arraySize, minElement, maxElement);


	ScalarEncoder encoder = ScalarEncoder(5, 10, 20, 20, true);

	UInt output_width = encoder.getOutputWidth();
	cout << "output_width:" << output_width << endl;

	vector<UInt> inputDimensions = {output_width}, columnDimensions = {50};

	SpatialPooler sp(inputDimensions, columnDimensions);
	TemporalMemory tm({ columnDimensions });

	UInt* inputVector = (UInt*)calloc(output_width, sizeof(UInt));
	UInt* columnVector = (UInt*)calloc(50, sizeof(UInt));

	for (int i = 0; i < arraySize; ++i) {
		cout << "iteration:" << i + 1 << endl;
		memset(inputVector, 0, sizeof(UInt)*output_width);
		encoder.encodeIntoArray(date[i], inputVector);
		printCodedDate(inputVector, output_width);
		sp.compute(inputVector, true, columnVector);
		printCodedDate(columnVector, 50);
		tm.compute(50, columnVector);
	}
	

}