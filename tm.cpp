#include <vector>
#define CONNECTED_PERMANENCE 1
#define ACTIVATION_THRESHOLD 1
#define LEARNING_THRESHOLD   1
#define LEARNING_ENABLED     1
#define PERMANENCE_INCREMENT 1
#define PERMANENCE_DECREMENT 1
#define SYNAPSE_SAMPLE_SIZE  1
#define PREDICTED_DECREMENT  1
#define INITIAL_PERMANENCE   1

using namespace std;

typedef int cell;

struct synapse
{
    int index_presynapticCell; // 保存细胞索引
    int permanence;
};

struct segment
{
    int index_cell; // 细胞索引
    vector<synapse> synapses;
    int numActivePotentialSynapses;
};



struct column
{
    vector<int> index_activeSegments;   // t-1时刻
    vector<int> index_matchingSegments; // t-1时刻
    vector<int> cells;
};



const int numColumns = 2048;

vector<column> columns;

vector<cell> cells; // t-1时刻, cells[i] == 1表示其状态为active
vector<int> cellsNumSegments; // t-1时刻和t时刻
vector<cell> activeCells; // t时刻
vector<cell> winnerCells; // t时刻

vector<segment> segments;
vector<segment> activeSegments;   // t时刻 TODO: 改成segments的索引节省内存
vector<segment> matchingSegments; // t时刻 TODO: 同上

void activateCells(vector<int> input, bool learn) {
    int sz = input.size();
    for (int i = 0; i < sz; ++i) {
        column &column = column[i];
        if (input[i] == true) {
            if (columns.index_activeSegments.size() > 0) activatePredictedColumn(column);
            else burstColumn(column);
        }
        else {
            if (columns.index_matchingSegments.size() > 0) punishPredictedColumn(column);
        }
    }
}

void activatePredictedColumn(column &column) {
    for (int index_segment : column.index_activeSegments) {
        segment &segment = segments[index_segment];
        activeCells(t).push_back(segment.index_cell);
        winnerCells(t).push_back(segment.index_cell);
        if (LEARNING_ENABLED) {
            for (synapse &synapse : segment.synapses) {
                if (cells(t-1)[synapse.index_presynapticCell] == 1)
                    synapse.permanence += PERMANENCE_INCREMENT;
                else
                    synapse.permanence -= PERMANENCE_DECREMENT;
            }
            int newSynapseCount = SYNAPSE_SAMPLE_SIZE - segment.numActivePotentialSynapses;
            growSynapses(segment, newSynapseCount); // TODO
        }
    }
}

void burstColumn(column &column) {
    for (cell cell : column.cells) {
        activeCells(t).push_back(cell);
        cells(t)[cell] = 1;
    }
    int learningSegment, winnerCell;
    if (column.index_matchingSegments.size() > 0) {
        learningSegment = bestMatchingSegment(column); // TODO
        winnerCell = learningSegment.cell;
    }
    else {
        winnerCell = leastUsedCell(column); // TODO
        if (LEARNING_ENABLED)
            learningSegment = growNewSegment(winnerCell); // TODO
    }
    winnerCells(t).push_back(winnerCell);
    if (LEARNING_ENABLED) {
        segment &segment = segments[learningSegment];
        for (synapse &synapse : segment.synapses) {
            if (cells(t-1)[synapse.index_presynapticCell] == 1)
                synapse.permanence += PERMANENCE_INCREMENT;
            else
                synapse.permanence -= PERMANENCE_DECREMENT;
            int newSynapseCount = SYNAPSE_SAMPLE_SIZE - segment.numActivePotentialSynapses;
            growSynapses(learningSegment, newSynapseCount);
        }
        
    }
}

void punishPredictedColumn(column &column) {
    if (LEARNING_ENABLED) {
        for (int index_segment : column.index_matchingSegments) {
            segment &segment = segments[index_segment];
            for (synapse &synapse : segment.synapses)
                if (cells(t-1)[synapse.index_presynapticCell] == 1)
                    synapse.permanence -= PREDICTED_DECREMENT;
        }
    }
}

void activateDendrites(bool learn) {
    for (segment &segment : segments) { // TODO: 修改为索引访问
        int numActiveConnected = 0;
        int numActivePotential = 0;
        for (synapse &synapse : segment.synapses) {
            if (cells(t)[synapse.index_presynapticCell]) {
                if (synapse.permanence >= CONNECTED_PERMANENCE) numActiveConnected++;
                if (synapse.permanence >= 0) numActivePotential++;
            }
        }
        
        if (numActiveConnected >= ACTIVATION_THRESHOLD) activeSegments(t).push_back(segment); // TODO
        if (numActivePotential >= LEARNING_THRESHOLD) matchingSegments(t).push_back(segment); // TODO
        segment.numActivePotentialSynapses = numActivePotential;
    }
}

cell& leastUsedCell(column &column) {
    int fewestSegments = INT_MAX;
    for (cell cell : column.cells)
        fewestSegments = min(fewestSegments, cellsNumSegments[cell]);
    vector<cell> leastUsedCells;
    for (cell cell : column.cells)
        if (cellsNumSegments[cell] == fewestSegments)
            leastUsedCells.push_back(cell);
    cell res = chooseRandom(leastUsedCells); // TODO
    return res;
}

segment& bestMatchingSegment(column &column) {
    segment res;
    int bestScore = -1;
    for (int index_matchingSegment : column.index_matchingSegments) {
        segment &segment = segments[index_matchingSegment];
        if (segment.numActivePotentialSynapses > bestScore) {
            res = segment;
            bestScore = segment.numActivePotentialSynapses;
        }
    }
    return res;
}

void growSynapses(segment &segment, int newSynapseCount) {
    vector<cell> candidates(winnerCells(t-1).begin(), winnerCells(t-1).end());
    while (candidates.size() > 0 && newSynapseCount > 0) {
        cell &presynapticCell = chooseRandom(candidates); // TODO
        candidates.remove(presynapticCell) // TODO
        bool alreadyConnected = false;
        for (synapse &synapse: segment.synapses)
            if (synapse.index_presynapticCell == presynapticCell)
                alreadyConnected = true;
        if (!alreadyConnected) {
             createNewSynapse(segment, presynapticCell, INITIAL_PERMANENCE);
             newSynapseCount--;
        }
    }
}