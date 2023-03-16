#include <random>
#include <algorithm>
#include <iostream>

#include "newsp.hpp"

using namespace std;

SpatialPooler::SpatialPooler(vector<UInt> inputDimensions, vector<UInt> columnDimensions,
                UInt potentialRadius, Real potentialPct,
                bool globalInhibition, Real localAreaDensity,
                UInt numActiveColumnsPerInhArea,
                UInt stimulusThreshold, Real synPermInactiveDec,
                Real synPermActiveInc, Real synPermConnected,
                Real minPctOverlapDutyCycles,
                UInt dutyCyclePeriod, Real boostStrength) {
    
    // 计算输入bit数
    inputDimensions_ = inputDimensions;
    numInputs_ = 1;
    for (auto &dim : inputDimensions) {
        numInputs_ *= dim;
    }

    // 计算columns数量
    columnDimensions_ = columnDimensions;
    numColumns_ = 1;
    for (auto &dim : columnDimensions) {
        numColumns_ *= dim;
    }

    potentialRadius_    = potentialRadius;
    potentialPct_       = potentialPct;
    globalInhibition_   = globalInhibition;
    localAreaDensity_   = localAreaDensity;
    stimulusThreshold_  = stimulusThreshold;
    synPermInactiveDec_ = synPermInactiveDec;
    synPermActiveInc_   = synPermActiveInc;
    synPermConnected_   = synPermConnected;
    dutyCyclePeriod_    = dutyCyclePeriod;
    boostStrength_      = boostStrength;
    minPctOverlapDutyCycles_    = minPctOverlapDutyCycles;
    numActiveColumnsPerInhArea_ = numActiveColumnsPerInhArea;

    // 为一些数组分配空间
    overlap_.assign(numColumns_, 0);
    boostedOverlaps_.assign(numColumns_, 0);
    boostFactors_.assign(numColumns_, 1.0);
    // neighbours_.assign(numColumns_, vector<UInt>(4));
    activeDutyCycle_.assign(numColumns_, 0);
    overlapDutyCycle_.assign(numColumns_, 0);
    minOverlapDutyCycles_.assign(numColumns_, 0);
    activeCount_.assign(numColumns_, 0);
    overlapCount_.assign(numColumns_, 0);

    /* 为每个column分配突触 */ 
    // assign receptive field
    mapIndex_(numColumns_, numInputs_, potentialRadius, potentialPools_);

    // assign synapses
    for (auto &field : potentialPools_) {
        assignSynapses_(field);
    }
    
    // compute inhibitionRadius
    updateInhibitionRadius_();

    // get neighbour
    updateNeighbour_();
}

void SpatialPooler::compute(UInt inputVector[], 
            bool learn, UInt activeVector[]) {
    // 更新dutyCycle时要用到
    iterationNum_++;
    // overlap
    cout << "overlap" << endl;
    calculateOverlap_(inputVector);
    if (learn) {
        cout << "boostOverlap_" << endl;

        boostOverlap_();
    }
    else {
        cout << "boostedOverlaps_" << endl;

        boostedOverlaps_.assign(overlap_.begin(), overlap_.end());
    }
    cout << "inhibitColumns_" << endl;

    inhibitColumns_(activeVector);


    if (learn) {
        cout << "adaptSynapses_" << endl;

        adaptSynapses_();
        cout << "updateDutyCycles_" << endl;

        updateDutyCycles_();
        cout << "updateBoostFactors_" << endl;

        updateBoostFactors_();
        cout << "bumpUpWeakColumns_" << endl;

        bumpUpWeakColumns_();
        if (iterationNum_ % dutyCyclePeriod_ == 0) {
            cout << "updateInhibitionRadius_" << endl;

            updateInhibitionRadius_();
            cout << "updateNeighbour_" << endl;

            updateNeighbour_();
            cout << "updateMinDutyCycles_" << endl;

            updateMinDutyCycles_();
        }
    }

}


void SpatialPooler::mapIndex_(const UInt sourceSize, const UInt targetSize, 
            UInt mapRadius, vector<vector<UInt>> &container) {
    container.clear();
    if (2 * mapRadius + 1 >= targetSize) {
        container.assign(sourceSize, {0, targetSize});
        return;
    }
    if ((mapRadius*2+1)*sourceSize < targetSize) {
        cout << "(mapRadius*2+1)*sourceSize < targetSize无法产生公平映射" << endl;
        throw logic_error("-1");
    }
    if (targetSize == sourceSize) {
        for (int i = 0; i < sourceSize; ++i) {
            vector<UInt> temp;
            int left = i - mapRadius, right = i + mapRadius;
            if (left < 0) {
                temp = {0, UInt(right) +1, sourceSize + left, sourceSize};
            }
            else if (right >= sourceSize) {
                temp = {UInt(left), sourceSize, 0, right-sourceSize+1};
            }
            else {
                temp = { UInt(left), UInt(right) + 1};
            }
            container.push_back(temp);
        }
    }
}


void SpatialPooler::assignSynapses_(vector<UInt> &receptiveField) {
    int numRange = receptiveField.size() / 2;
    vector<int> potentialField;
    for (int i = 0; i < numRange; ++i) {
        for (int j = receptiveField[i*2]; j < receptiveField[i*2+1]; ++j) {
            potentialField.push_back(j);
        }
    }
    // randomly get numPotential 
    // 以下算法无法确定正确性
    vector<synapse> potentialSyn, connectedSyn;
    for (auto &idx : potentialField) {
        if (rand() % 10 / 10.0 > 0.5) continue;
        synapse syn;
        syn.sourceInputIndex = idx;
        syn.permanence = rand() % 11 / 10.0;
        if (syn.permanence >= synPermConnected_) {
            connectedSyn.push_back(syn);
        }
        potentialSyn.push_back(syn);
    }
    potentialSyn_.push_back(potentialSyn);
    connectedSyn_.push_back(connectedSyn);
}

void SpatialPooler::updateInhibitionRadius_() {
    if (globalInhibition_) {
        inhibitionRadius_ = numColumns_;
        return;
    }
    UInt numConnectedSyn = 0;
    for (auto &col : connectedSyn_) {
        numConnectedSyn += col.size();
    }
    inhibitionRadius_ = numConnectedSyn / numColumns_;
    inhibitionRadius_ = inhibitionRadius_ > numColumns_ ? \
                            numColumns_ : inhibitionRadius_;
}

void SpatialPooler::calculateOverlap_(UInt inputVector[]) {
    for (int i = 0; i < numColumns_; ++i) {
        int tmp = 0;
        for (auto &syn : connectedSyn_[i]) {
            tmp += inputVector[syn.sourceInputIndex];
        }
        overlap_[i] = tmp;
    }
}

void SpatialPooler::boostOverlap_() {
    for (int i = 0; i < numColumns_; ++i) {
        boostedOverlaps_[i] = overlap_[i] * boostFactors_[i];
    }
}

void SpatialPooler::updateNeighbour_() {
    mapIndex_(numColumns_, numColumns_, inhibitionRadius_, neighbours_);
}

void SpatialPooler::inhibitColumns_(UInt activeVector[]) {
    int numActive = localAreaDensity_ > 0 ? (inhibitionRadius_*2+1)
         * localAreaDensity_ : numActiveColumnsPerInhArea_;
    activeColumns_.clear();
    Real globalNumActive = kthScore(0, numActive);
    for (int i = 0; i < numColumns_; ++i) {
        Real minLocalActivity = 
            globalInhibition_ == true ? globalNumActive : kthScore(i, numActive);
        if (boostedOverlaps_[i] > stimulusThreshold_) {
            overlapCount_[i]++;
            if (boostedOverlaps_[i] >= minLocalActivity) {
                activeColumns_.push_back(i);
                activeCount_[i]++;
                activeVector[i] = 1;
            }
        }
        else activeVector[i] = 0;
    }
}

Real SpatialPooler::kthScore(UInt colIndex, UInt numActive) {
    vector<Real> neighboursOverlap;
    int numRange = neighbours_[colIndex].size() / 2;
    for (int i = 0; i < numRange; ++i) {
        for (int j = neighbours_[colIndex][i*2]; j < neighbours_[colIndex][i*2+1]; ++j) {
            neighboursOverlap.push_back(boostedOverlaps_[j]);
        }
    }
    sort(neighboursOverlap.begin(), neighboursOverlap.end());
    return neighboursOverlap[numActive-1];
}

void SpatialPooler::adaptSynapses_() {
    for (auto columnIndex : activeColumns_) {
        for (auto &syn : potentialSyn_[columnIndex]) {
            if (syn.permanence > synPermConnected_) {
                syn.permanence += synPermActiveInc_;
                syn.permanence = min(Real(1), syn.permanence);
            }
            else {
                syn.permanence -= synPermInactiveDec_;
                syn.permanence = max(Real(0), syn.permanence);
            }
        }
    }
}

void SpatialPooler::updateDutyCycles_() {
    for (int i = 0; i < numColumns_; ++i) {
        activeDutyCycle_[i] = activeCount_[i] / iterationNum_;
        overlapDutyCycle_[i] = overlapCount_[i] / iterationNum_;
    }
}

void SpatialPooler::updateBoostFactors_() {
    Real activeDutyCycleNeighbors = 0;
    for (int idx = 0; idx < numColumns_; ++idx) {
        int numRange = neighbours_[idx].size() / 2;
        Real sumDutyCycle = 0;
        for (int i = 0; i < numRange; ++i) {
            for (int j = neighbours_[idx][i*2]; j < neighbours_[idx][i*2+1]; ++j) {
                sumDutyCycle += activeDutyCycle_[j];
            }
        }
        activeDutyCycleNeighbors = sumDutyCycle / (2* inhibitionRadius_+1);
        boostFactors_[idx] = 
            exp(boostStrength_ * (activeDutyCycleNeighbors - activeDutyCycle_[idx]));
    }
    
}

void SpatialPooler::bumpUpWeakColumns_() {
    for (int i = 0; i < numColumns_; ++i) {
        if (overlapDutyCycle_[i] < minOverlapDutyCycles_[i]) {
            increasePermanences_(i, 0.1);
        }
    }
}

void SpatialPooler::updateMinDutyCycles_() {
    if (globalInhibition_) {
        Real maxOverlapDutyCycles =
            *max_element(overlapDutyCycle_.begin(), overlapDutyCycle_.end());
        fill(minOverlapDutyCycles_.begin(), minOverlapDutyCycles_.end(),
            minPctOverlapDutyCycles_ * maxOverlapDutyCycles);
        return;
    }
    // 非全局抑制更新MinDutyCycles，效率很低
    for (int i = 0; i < numColumns_; ++i) {
        vector<Real> neighbourOverlapDutyCycles;
        int numRange = neighbours_[i].size() / 2;
        for (int i = 0; i < numRange; ++i) {
            for (int j = neighbours_[i][i*2]; j < neighbours_[i][i*2+1]; ++j) {
                neighbourOverlapDutyCycles.push_back(overlapDutyCycle_[j]);
            }
        }
        Real maxOverlapDutyCycles =
            *max_element(neighbourOverlapDutyCycles.begin(), neighbourOverlapDutyCycles.end());
        minOverlapDutyCycles_[i] = minPctOverlapDutyCycles_ * maxOverlapDutyCycles;
    }
}

void SpatialPooler::increasePermanences_(UInt columnIndex, Real incFactor) {
    for (auto &syn : potentialSyn_[columnIndex]) {
        syn.permanence += synPermConnected_ * incFactor;
    }
}


