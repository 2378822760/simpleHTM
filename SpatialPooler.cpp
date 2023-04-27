// #define DEBUG
#include <random>
#include <algorithm>
#include <iostream>

#include "SpatialPooler.hpp"

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
#ifdef DEBUG
    cout << "compute" << endl;
#endif // DEBUG
    // 更新dutyCycle时要用到
    iterationNum_++;

    calculateOverlap_(inputVector);

    if (learn) boostOverlap_();
    else boostedOverlaps_.assign(overlap_.begin(), overlap_.end());
    inhibitColumns_(activeVector);

    if (learn) {
        adaptSynapses_();
        updateDutyCycles_();
        updateBoostFactors_();
        bumpUpWeakColumns_();
        if (iterationNum_ % dutyCyclePeriod_ == 0) {
            updateInhibitionRadius_();
            updateNeighbour_();
            updateMinDutyCycles_();
        }
    }

}

void SpatialPooler::mapIndexE_(const UInt beg, const UInt end, const UInt size,
            const UInt mapRadius, vector<vector<UInt>>& container) {
#ifdef DEBUG
    cout << "mapIndexE_" << endl;
#endif // DEBUG
    for (int i = beg; i < end; ++i) {
        vector<UInt> temp;
        int left = i - mapRadius, right = i + mapRadius;
        if (left < 0) {
            temp = { 0, UInt(right) + 1, size + left, size };
        }
        else if (right >= size) {
            temp = { UInt(left), size, 0, right - size + 1 };
        }
        else {
            temp = { UInt(left), UInt(right) + 1 };
        }
        container.push_back(temp);
    }
}

void SpatialPooler::mapIndex_(const UInt sourceSize, const UInt targetSize, 
            UInt mapRadius, vector<vector<UInt>> &container) {
#ifdef DEBUG
    cout << "mapIndex_" << endl;
#endif // DEBUG
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
        mapIndexE_(0, targetSize, targetSize, mapRadius, container);
    }
    else if (targetSize < sourceSize) {
        int coverCount = sourceSize / targetSize;
        int remainCount = sourceSize % targetSize;
        vector<vector<UInt>> temp;
        mapIndexE_(0, targetSize, targetSize, mapRadius, temp);
        for (int i = 0; i < coverCount; ++i) {
            container.insert(container.end(), temp.begin(), temp.end());
        }
        UInt beg = (targetSize - remainCount) / 2;
        mapIndexE_(beg, beg + remainCount, targetSize, mapRadius, container);
    }
}


void SpatialPooler::assignSynapses_(vector<UInt> &receptiveField) {
#ifdef DEBUG
    cout << "assignSynapses_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "updateInhibitionRadius_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "calculateOverlap_" << endl;
#endif // DEBUG
    for (int i = 0; i < numColumns_; ++i) {
        int tmp = 0;
        for (auto &syn : connectedSyn_[i]) {
            tmp += inputVector[syn.sourceInputIndex];
        }
        overlap_[i] = tmp;
    }
}

void SpatialPooler::boostOverlap_() {
#ifdef DEBUG
    cout << "boostOverlap_" << endl;
#endif // DEBUG
    for (int i = 0; i < numColumns_; ++i) {
        boostedOverlaps_[i] = overlap_[i] * boostFactors_[i];
    }
}

void SpatialPooler::updateNeighbour_() {
#ifdef DEBUG
    cout << "updateNeighbour_" << endl;
#endif // DEBUG
    mapIndex_(numColumns_, numColumns_, inhibitionRadius_, neighbours_);
}

void SpatialPooler::inhibitColumns_(UInt activeVector[]) {
#ifdef DEBUG
    cout << "inhibitColumns_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "kthScore" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "adaptSynapses_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "updateDutyCycles_" << endl;
#endif // DEBUG
    for (int i = 0; i < numColumns_; ++i) {
        activeDutyCycle_[i] = activeCount_[i] / iterationNum_;
        overlapDutyCycle_[i] = overlapCount_[i] / iterationNum_;
    }
}

void SpatialPooler::updateBoostFactors_() {
#ifdef DEBUG
    cout << "updateBoostFactors_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "bumpUpWeakColumns_" << endl;
#endif // DEBUG
    for (int i = 0; i < numColumns_; ++i) {
        if (overlapDutyCycle_[i] < minOverlapDutyCycles_[i]) {
            increasePermanences_(i, 0.1);
        }
    }
}

void SpatialPooler::updateMinDutyCycles_() {
#ifdef DEBUG
    cout << "updateMinDutyCycles_" << endl;
#endif // DEBUG
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
#ifdef DEBUG
    cout << "increasePermanences_" << endl;
#endif // DEBUG
    for (auto &syn : potentialSyn_[columnIndex]) {
        syn.permanence += synPermConnected_ * incFactor;
    }
}


