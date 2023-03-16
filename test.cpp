#include "encoder.hpp"
#include "newsp.hpp"

#include <iostream>
#include <cstdio>



using namespace std;

/**
 * ����arr�����½磬�Ͻ籣����max���½籣����min
 */
void getArrayRange(const double arr[], int arraySize, double& min, double& max) {
	// ��ƺ���ͬʱ�����С���ֵ
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
	double date[] = { 10,10.5,11,11.5,12,12.5,13,14,15,16,17,18,19 };
	int arraySize = sizeof(date) / sizeof(double);
	double minElement, maxElement;

	getArrayRange(date, arraySize, minElement, maxElement);

	// ����ʵ��
	ScalarEncoder encoder = ScalarEncoder(5, 10, 20, 20, true);
	// ��ȡ�������
	UInt output_width = encoder.getOutputWidth();
	cout << "output_width:" << output_width << endl;

	// ��ʼ��SP
	vector<UInt> inputDimensions = {output_width}, columnDimensions = {output_width};

	SpatialPooler sp(inputDimensions, columnDimensions);
	// �洢���ݱ���������
	UInt* inputVector = (UInt*)calloc(output_width, sizeof(UInt));
	UInt* columnVector = (UInt*)calloc(output_width, sizeof(UInt));
	// ��ʼ��һ���ռ��
	for (int i = 0; i < arraySize; ++i) {
		cout << "iteration:" << i + 1 << endl;
		memset(inputVector, 0, sizeof(UInt)*output_width);
		encoder.encodeIntoArray(date[i], inputVector);
		printCodedDate(inputVector, output_width);
		sp.compute(inputVector, true, columnVector);
		printCodedDate(columnVector, output_width);
	}
	

}