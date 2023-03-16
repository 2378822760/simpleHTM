#include "encoder.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

ScalarEncoder::ScalarEncoder(int w, double minValue, double maxValue, int bucketNum, bool clipInput)
	: w_(w), minValue_(minValue), maxValue_(maxValue), bucketNum_(bucketNum), clipInput_(clipInput) {

	if (w <= 0) {
		cout << "w must be > 0" << endl;
		throw invalid_argument("-1");
	}

	if (bucketNum <= 0) {
		cout << "bucketNum = " << bucketNum << "must > 0" << endl;
		throw invalid_argument("-1");
	}

	if (minValue >= maxValue) {
		cout << "minValue must be < maxValue. minValue=" << minValue
			<< " maxValue=" << maxValue << endl;
		throw invalid_argument("-1");
	}

	// ¼ÆËãÊä³ö¿í¶È
	n_ = bucketNum + w - 1;

}

void ScalarEncoder::encodeIntoArray(double input, UInt output[]) {
	if (input < minValue_ || input > maxValue_) {
		if (clipInput_) {
			input = input < minValue_ ? minValue_ : maxValue_;
		}
		else {
			cout << "input (" << input << ") out of range [" << minValue_
				<< ", " << maxValue_ << "]";
			throw out_of_range("-1");
		}
	}

	const int iBucket = round(bucketNum_ * (input - minValue_) / (maxValue_ - minValue_));
	// cout << iBucket << endl;
	int endedIdx = w_ + iBucket;

	fill(output + iBucket, output + endedIdx, 1);

}
