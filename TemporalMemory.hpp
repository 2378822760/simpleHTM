//
// Created by 曹 on 2023/4/24.
//

#ifndef TEMPORALMEMORY_HPP
#define TEMPORALMEMORY_HPP

#endif // TEMPORALMEMORY_HPP

#include "type.hpp"
#include <vector>
#include <unordered_set>

using namespace std;

struct Synapse {
    CellIndex presynapticCell;
    Permanence permanence;
};

struct Segment {
    CellIndex cell;             // 保存该seg所依附的细胞
    vector<Synapse> synapses;   // syn不能单独存在，要依附于seg
    UInt numActivePotentialSynapses;
};

struct Column {
    vector<SegIndex> activeSegments;   // t-1时刻，每当一次计算后由activeSegments更新
    vector<SegIndex> matchingSegments; // t-1时刻，每当一次计算后由matchingSegments更新
    vector<CellIndex> cells;
};

class TemporalMemory {
private:
    /* 全局变量 */
    UInt ACTIVATION_THRESHOLD{};
    UInt LEARNING_THRESHOLD{};
    UInt SYNAPSE_SAMPLE_SIZE{};
    Permanence PERMANENCE_INCREMENT{};
    Permanence PERMANENCE_DECREMENT{};
    Permanence CONNECTED_PERMANENCE{};
    Permanence PREDICTED_DECREMENT{};
    Permanence INITIAL_PERMANENCE{};

    /* 过程变量 */
    bool LEARNING_ENABLED{};
    UInt iteration_{};

    /* 初始化变量，在构造函数或者调用initialize函数中初始化 */
    vector<UInt> columnDimensions_;
    UInt cellsPerColumn_{};
    UInt numColumns_{};

    vector<Column> columns;
    vector<ColumnIndex> cellsFa;    // 纽带，保存每个cell所依附的col
    unordered_set<CellIndex> candidates;   // t-1时刻，每当一次计算后由winnerCells更新
    vector<bool> cells;             // t-1时刻，每当一次计算后由activeCells更新
    vector<bool> winnerCells;       // t时刻
    vector<bool> activeCells;       // t时刻
    vector<UInt> cellsNumSegments;  // 全局共享，无需区分时刻

    vector<Segment> segments;
    vector<SegIndex> activeSegments;   // t时刻
    vector<SegIndex> matchingSegments; // t时刻
public:
    TemporalMemory();

    explicit TemporalMemory(const vector<UInt>& columnDimensions, UInt cellsPerColumn = 32,
        UInt activationThreshold = 13, Permanence initialPermanence = 0.21,
        Permanence connectedPermanence = 0.50, UInt LearningThreshold = 10,
        Permanence permanenceIncrement = 0.10, Permanence permanenceDecrement = 0.10,
        Permanence predictedSegmentDecrement = 0.0);

    ~TemporalMemory();

    void initialize(const vector<UInt>& columnDimensions, UInt cellsPerColumn = 32,
        UInt activationThreshold = 13, Permanence initialPermanence = 0.21,
        Permanence connectedPermanence = 0.50, UInt LearningThreshold = 10,
        Permanence permanenceIncrement = 0.10, Permanence permanenceDecrement = 0.10,
        Permanence predictedSegmentDecrement = 0.0);

    void compute(size_t activeColumnsSize, const UInt activeColumns[],
        bool learn = true);

    void activateCells(size_t activeColumnsSize, const UInt activeColumns[]);

    void activatePredictedColumn(Column& column);

    void burstColumn(Column& column);

    void punishPredictedColumn(Column& column);

    void activateDendrites();

    /* Helper Function */

    SegIndex growNewSegment(CellIndex cell);

    CellIndex leastUsedCell(Column& column);

    SegIndex bestMatchingSegment(Column& column);

    void growSynapses(SegIndex segment, UInt newSynapseCount);

    template<class T>
    T chooseRandom(vector<T>& Candidates);

    template<class T>
    T chooseRandom(unordered_set<T>& Candidates);

    void createNewSynapse(SegIndex index_seg, CellIndex presynapticCell);

    void updateSegments();

    void updateCells();
};