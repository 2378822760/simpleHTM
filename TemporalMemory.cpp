//
// Created by 曹 on 2023/4/25.
//
// #define DEBUG
#include <climits>
#include <ctime>
#include <cstdlib>
#include <numeric>
#include <iostream>
#include "TemporalMemory.hpp"

TemporalMemory::TemporalMemory() = default;

TemporalMemory::TemporalMemory(const vector<UInt>& columnDimensions, UInt cellsPerColumn, UInt activationThreshold,
    Permanence initialPermanence, Permanence connectedPermanence, UInt LearningThreshold,
    Permanence permanenceIncrement, Permanence permanenceDecrement,
    Permanence predictedSegmentDecrement) {
    initialize(columnDimensions, cellsPerColumn, activationThreshold,
        initialPermanence, connectedPermanence, LearningThreshold,
        permanenceIncrement, permanenceDecrement,
        predictedSegmentDecrement);
}

TemporalMemory::~TemporalMemory() = default;

void TemporalMemory::initialize(const vector<UInt>& columnDimensions, UInt cellsPerColumn, UInt activationThreshold,
    Permanence initialPermanence, Permanence connectedPermanence, UInt LearningThreshold,
    Permanence permanenceIncrement, Permanence permanenceDecrement,
    Permanence predictedSegmentDecrement) {

#ifdef DEBUG
    cout << "initialize" << "\n";
#endif // DEBUG

    /* TODO: 编写参数合法性检测 */

    numColumns_ = 1;
    for (UInt n : columnDimensions) {
        columnDimensions_.push_back(n);
        numColumns_ *= n;
    }

    cellsPerColumn_ = cellsPerColumn;
    ACTIVATION_THRESHOLD = activationThreshold;
    INITIAL_PERMANENCE = initialPermanence;
    CONNECTED_PERMANENCE = connectedPermanence;
    LEARNING_THRESHOLD = LearningThreshold;
    PERMANENCE_INCREMENT = permanenceIncrement;
    PERMANENCE_DECREMENT = permanenceDecrement;
    PREDICTED_DECREMENT = predictedSegmentDecrement;

    iteration_ = 0;


    columns.clear();
    cells.clear();
    cellsNumSegments.clear();
    segments.clear();
    activeCells.clear();
    winnerCells.clear();
    activeSegments.clear();
    matchingSegments.clear();
    columns.resize(numColumns_);
    cells.resize(numColumns_ * cellsPerColumn_);
    cellsNumSegments.resize(numColumns_ * cellsPerColumn_);
    activeCells.resize(numColumns_ * cellsPerColumn_);
    winnerCells.resize(numColumns_ * cellsPerColumn_);
    cellsFa.resize(numColumns_ * cellsPerColumn_);
    for (int i = 0; i < numColumns_; ++i) { // 声明每个col的所有cell
        Column& column = columns[i];
        column.cells.resize(cellsPerColumn_);
        iota(column.cells.begin(), column.cells.end(), i * cellsPerColumn_);
        auto beg = cellsFa.begin() + i * cellsPerColumn_;
        auto end = beg + cellsPerColumn_;
        fill(beg, end, i);
    }

}

void TemporalMemory::compute(std::size_t activeColumnsSize, const UInt* activeColumns, bool learn) {
    LEARNING_ENABLED = learn;
    activateCells(activeColumnsSize, activeColumns);
    activateDendrites();
    updateSegments();
    updateCells();
}

void TemporalMemory::activateCells(std::size_t activeColumnsSize, const UInt* activeColumns) {
#ifdef DEBUG
    cout << "activateCells" << "\n";
#endif // DEBUG
    for (int i = 0; i < activeColumnsSize; ++i) {
        Column& column = columns[i];
        if (activeColumns[i] == true) {
            if (column.activeSegments.size() > 0) activatePredictedColumn(column);
            else burstColumn(column);
        }
        else {
            if (column.activeSegments.size() > 0) punishPredictedColumn(column);
        }
    }
}

void TemporalMemory::activatePredictedColumn(Column& column) {
#ifdef DEBUG
    cout << "activatePredictedColumn" << "\n";
#endif // DEBUG
    for (SegIndex n : column.activeSegments) {
        Segment& segment = segments[n];
        activeCells[segment.cell] = true; // 新一轮循环之前应保证该列表为空
        winnerCells[segment.cell] = true; // t时刻
        if (LEARNING_ENABLED) {
            for (Synapse& synapse : segment.synapses) {
                if (cells[synapse.presynapticCell] == 1) // 直到该操作前cells表不应被更新 @time t-1时刻
                    synapse.permanence += PERMANENCE_INCREMENT;
                else
                    synapse.permanence -= PERMANENCE_DECREMENT;
            }
            UInt newSynapseCount = SYNAPSE_SAMPLE_SIZE - segment.numActivePotentialSynapses;
            growSynapses(n, newSynapseCount);
        }
    }
}

void TemporalMemory::burstColumn(Column& column) {
#ifdef DEBUG
    cout << "burstColumn" << "\n";
#endif // DEBUG
    for (CellIndex cell : column.cells) {
        activeCells[cell] = true;
    }
    SegIndex learningSegment;
    CellIndex winnerCell;
    if (column.matchingSegments.size() > 0) {
        learningSegment = bestMatchingSegment(column);
        winnerCell = segments[learningSegment].cell;
    }
    else {
        winnerCell = leastUsedCell(column);
        if (LEARNING_ENABLED)
            learningSegment = growNewSegment(winnerCell);
    }
    winnerCells[winnerCell] = true;
    if (LEARNING_ENABLED) {
        Segment& segment = segments[learningSegment];
        for (Synapse& synapse : segment.synapses) {
            if (cells[synapse.presynapticCell] == 1) // cells(t-1)在此步骤不应该被重置 @time t-1
                synapse.permanence += PERMANENCE_INCREMENT;
            else
                synapse.permanence -= PERMANENCE_DECREMENT;
            UInt newSynapseCount = SYNAPSE_SAMPLE_SIZE - segment.numActivePotentialSynapses;
            growSynapses(learningSegment, newSynapseCount);
        }
    }
}

SegIndex TemporalMemory::growNewSegment(CellIndex cell) {
#ifdef DEBUG
    cout << "growNewSegment CellIndex: " << cell << \
        "cellsNumSegments.size: " << cellsNumSegments.size() << "\n";
#endif // DEBUG
    SegIndex newSegment = segments.size();
    segments.push_back({ cell });
    cellsNumSegments[cell]++;
    return newSegment;
}

void TemporalMemory::punishPredictedColumn(Column& column) {
#ifdef DEBUG
    cout << "punishPredictedColumn" << "\n";
#endif // DEBUG
    if (LEARNING_ENABLED) {
        for (SegIndex index_segment : column.matchingSegments) {
            Segment& segment = segments[index_segment];
            for (Synapse& synapse : segment.synapses)
                if (cells[synapse.presynapticCell] == 1) // @time t-1
                    synapse.permanence -= PREDICTED_DECREMENT;
        }
    }
}

void TemporalMemory::activateDendrites() {
#ifdef DEBUG
    cout << "activateDendrites" << "\n";
#endif // DEBUG
    for (int i = 0; i < segments.size(); ++i) {
        Segment& segment = segments[i];
        UInt numActiveConnected = 0;
        UInt numActivePotential = 0;
        for (Synapse& synapse : segment.synapses) {
            if (activeCells[synapse.presynapticCell]) { // @time t
                if (synapse.permanence >= CONNECTED_PERMANENCE) numActiveConnected++;
                if (synapse.permanence >= 0) numActivePotential++;
            }
        }
        if (numActiveConnected >= ACTIVATION_THRESHOLD) activeSegments.push_back(i); // @time t
        if (numActivePotential >= LEARNING_THRESHOLD) matchingSegments.push_back(i); // @time t
        segment.numActivePotentialSynapses = numActivePotential;
    }
}

CellIndex TemporalMemory::leastUsedCell(Column& column) {
#ifdef DEBUG
    cout << "leastUsedCell" << "\n";
#endif // DEBUG
    UInt fewestSegments = INT_MAX;
    for (CellIndex cell : column.cells)
        fewestSegments = min(fewestSegments, cellsNumSegments[cell]);
    vector<CellIndex> leastUsedCells;
    for (CellIndex cell : column.cells)
        if (cellsNumSegments[cell] == fewestSegments)
            leastUsedCells.push_back(cell);
    CellIndex res = chooseRandom(leastUsedCells);
    return res;
}

SegIndex TemporalMemory::bestMatchingSegment(Column& column) {
#ifdef DEBUG
    cout << "bestMatchingSegment" << "\n";
#endif // DEBUG
    SegIndex res;
    int bestScore = -1;
    for (SegIndex index_matchingSegment : column.matchingSegments) {
        Segment& segment = segments[index_matchingSegment];
        if (segment.numActivePotentialSynapses > bestScore) {
            res = index_matchingSegment;
            bestScore = segment.numActivePotentialSynapses;
        }
    }
    return res;
}

void TemporalMemory::growSynapses(SegIndex index_segment, UInt newSynapseCount) {
#ifdef DEBUG
    cout << "growSynapses" << "\n";
#endif // DEBUG
    Segment& segment = segments[index_segment];
    while (candidates.size() > 0 && newSynapseCount > 0) {
        CellIndex presynapticCell = chooseRandom(candidates);
        candidates.erase(presynapticCell); // TODO: 如何快速从vector中删除已用的candidate
        bool alreadyConnected = false;
        for (Synapse& synapse : segment.synapses)
            if (synapse.presynapticCell == presynapticCell)
                alreadyConnected = true;
        if (!alreadyConnected) {
            createNewSynapse(index_segment, presynapticCell);
            newSynapseCount--;
        }
    }
}

void TemporalMemory::createNewSynapse(SegIndex index_seg, CellIndex presynapticCell) {
#ifdef DEBUG
    cout << "createNewSynapse" << "\n";
#endif // DEBUG
    Segment& segment = segments[index_seg];
    segment.synapses.push_back({ presynapticCell, INITIAL_PERMANENCE });
}

void TemporalMemory::updateSegments() {
#ifdef DEBUG
    cout << "updateSegments" << "\n";
#endif // DEBUG
    for (auto& column : columns) {
        column.matchingSegments.clear();
        column.activeSegments.clear();
    }
    for (auto index_seg : matchingSegments)
        columns[cellsFa[segments[index_seg].cell]].matchingSegments.push_back(index_seg);

    for (auto index_seg : activeSegments)
        columns[cellsFa[segments[index_seg].cell]].activeSegments.push_back(index_seg);

}

void TemporalMemory::updateCells() {
#ifdef DEBUG
    cout << "updateCells" << "\n";
#endif // DEBUG
    int beg = int(cells.size()) - 1;
    swap(cells, activeCells);
    fill(activeCells.begin(), activeCells.end(), 0);
    candidates.clear();
    for (int i = beg; i >= 0; i--)
        if (winnerCells[i]) candidates.insert(i);
    fill(winnerCells.begin(), winnerCells.end(), 0);
}

template<class T>
T TemporalMemory::chooseRandom(unordered_set<T>& Candidates) {
#ifdef DEBUG
    cout << "chooseRandom set" << "\n";
#endif // DEBUG
    UInt seed = time(nullptr);
    srand(seed);
    size_t idx = ::rand() % Candidates.size();
    auto it = Candidates.begin();
    while (it != Candidates.end() && idx > 0) it++;
    return *it;
}

template<class T>
T TemporalMemory::chooseRandom(vector<T>& Candidates) {
#ifdef DEBUG
    cout << "chooseRandom vector" << "\n";
#endif // DEBUG
    UInt seed = time(nullptr);
    srand(seed);
    return Candidates[::rand() % Candidates.size()];
}