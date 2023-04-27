//
// Created by �� on 2023/4/24.
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
    CellIndex cell;             // �����seg��������ϸ��
    vector<Synapse> synapses;   // syn���ܵ������ڣ�Ҫ������seg
    UInt numActivePotentialSynapses;
};

struct Column {
    vector<SegIndex> activeSegments;   // t-1ʱ�̣�ÿ��һ�μ������activeSegments����
    vector<SegIndex> matchingSegments; // t-1ʱ�̣�ÿ��һ�μ������matchingSegments����
    vector<CellIndex> cells;
};

class TemporalMemory {
private:
    /* ȫ�ֱ��� */
    UInt ACTIVATION_THRESHOLD{};
    UInt LEARNING_THRESHOLD{};
    UInt SYNAPSE_SAMPLE_SIZE{};
    Permanence PERMANENCE_INCREMENT{};
    Permanence PERMANENCE_DECREMENT{};
    Permanence CONNECTED_PERMANENCE{};
    Permanence PREDICTED_DECREMENT{};
    Permanence INITIAL_PERMANENCE{};

    /* ���̱��� */
    bool LEARNING_ENABLED{};
    UInt iteration_{};

    /* ��ʼ���������ڹ��캯�����ߵ���initialize�����г�ʼ�� */
    vector<UInt> columnDimensions_;
    UInt cellsPerColumn_{};
    UInt numColumns_{};

    vector<Column> columns;
    vector<ColumnIndex> cellsFa;    // Ŧ��������ÿ��cell��������col
    unordered_set<CellIndex> candidates;   // t-1ʱ�̣�ÿ��һ�μ������winnerCells����
    vector<bool> cells;             // t-1ʱ�̣�ÿ��һ�μ������activeCells����
    vector<bool> winnerCells;       // tʱ��
    vector<bool> activeCells;       // tʱ��
    vector<UInt> cellsNumSegments;  // ȫ�ֹ�����������ʱ��

    vector<Segment> segments;
    vector<SegIndex> activeSegments;   // tʱ��
    vector<SegIndex> matchingSegments; // tʱ��
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