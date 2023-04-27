#include <vector>

#include "type.hpp"
/**
 * 突触包含了永久值permanence和连接的输入空间索引
 */
typedef struct synapse
{
	Real permanence;
	UInt sourceInputIndex;
}synapse;

using namespace std;


class SpatialPooler {
public:

    /**
     * 空间池的构造函数
     * @param inputDimensions 输入维度。用于计算输入的比特数
     * @param columnDimensions column维度。用于计算空间池column数量
     * @param potentialRadius 每个column的感受野半径，感受野大小为
     *      2 * potentialRadius + 1
     * @param potentialPct 每个column可以连接感受野内输入的比例。通过计算
     *      我们的得到一个column的隐含突触数量
     *      (2 * potentialRadius + 1) * potentialPct
     * @param globalInhibition 是否开启全局抑制。全局抑制意味着每个column
     *      的竞争对象为所有column。
     * @param localAreaDensity 在一个抑制半径内可以获胜的column的比例。该
     *      参数和numActiveColumnsPerInhArea只可以选择一个。若要使用该参数
     *      numActiveColumnsPerInhArea应为负数。
     * @param numActiveColumnsPerInhArea 在一个抑制半径内可以获胜的column
     *      的数量。该参数和localAreaDensity只可以选择一个。若要使用该参数
     *      localAreaDensity应为负数。
     * @param stimulusThreshold 刺激阈值。参与竞争的column的overlap值大于
     *      该阈值才可能被激活。
     * @param synPermInactiveDec 针对被激活的column，若其突触未连接，那么
     *      该突触的永久值 - synPermInactiveDec
     * @param synPermActiveInc 针对被激活的column，若其突触连接，那么
     *      该突触的永久值 + synPermActiveInc
     * @param synPermConnected 突触激活阈值。当突触的永久值大于该阈值时，
     *      该突触被激活。
     * @param minPctOverlapDutyCycles 用于更新minOverlapDutyCycles。
     * @param dutyCyclePeriod 更新dutyCycles时的迭代次数收敛于该参数。换而
     *      言之，如果当前迭代次数少于该值，那么计算dutyCycles时就除以迭代
     *      次数。当迭代次数大于该参数时，除以该参数。
     * @param boostStrength 更新boostFactors时，boostFunction使用的参数β。
     */
    SpatialPooler(vector<UInt> inputDimensions, vector<UInt> columnDimensions,
                UInt potentialRadius = 3, Real potentialPct = 0.5,
                bool globalInhibition = true, Real localAreaDensity = -1.0,
                UInt numActiveColumnsPerInhArea = 10,
                UInt stimulusThreshold = 0, Real synPermInactiveDec = 0.008,
                Real synPermActiveInc = 0.05, Real synPermConnected = 0.1,
                Real minPctOverlapDutyCycles = 0.001,
                UInt dutyCyclePeriod = 1000, Real boostStrength = 0.0);
    
    ~SpatialPooler() {};

    /**
     * 空间池最主要的函数。根据传递给函数的输入，计算出对应的激活column。开启
     *      学习后会进行参数更新。
     * @param inputVector 二值(0,1)向量。该向量通过编码器获得。向量长度应严
     *      格等于numInputs。
     * @param learn 如果开启学习，那么SP中一些参数会进行更新。
     * @param activeVector 二值(0,1)向量，若该值为1，那么说明该column被激活。
     *      该向量作为本次输入的输出。其大小为numColumns。
     */
    void compute(UInt inputVector[], bool learn, UInt activeVector[]);

    /**
     * 为每一个源索引映射一组目标索引（mapRadius为映射半径）。将每组索引存
     * 放到二维数组container。存放之前会先对该容器进行clear操作。
     * 其中每个一维向量的内容为{beg1, end1,(beg2, end2)}。最多包含两组索引
     * 范围，有第二组是因为边界元素索引为中心所以会超过左边界从右边界索引。
     */
    void mapIndex_(const UInt sourceSize, const UInt targetSize, UInt mapRadius,
        vector<vector<UInt>> &container);

    void mapIndexE_(const UInt beg, const UInt end, const UInt size,
        const UInt mapRadius, vector<vector<UInt>>& container);
    /**
     * （初始化时）为每个column分配突触。给定一个感受野范围（由mapIndex获得）
     * 根据参数potentialPct为感受野内输入随机分配突触，并随机赋初始值。
     * @param receptiveField 一维向量的内容为{beg1, end1,(beg2, end2)}，
     *      最多包含两组索引范围。
    */
    void assignSynapses_(vector<UInt> &receptiveField);

    /**
     * 更新半径。如果开启全局抑制，那么抑制范围就为numColumns。否则抑制半径
     * 需要更新为每个column连接突触个数的平均值。
    */
    void updateInhibitionRadius_();

    /**
     * 为每个输入计算overlap。每个column的overlap值为connectedSyn连接到为
     * 1的输入的数量。
     * @param inputVector 输入向量(二值0，1)，通过encoder获得。
    */
    void calculateOverlap_(UInt inputVector[]);

    /**
     * 开启学习后，每个column的overlap乘一个boostFactor，来作为最终的overlap。
    */
    void boostOverlap_();

    /**
     * 抑制。在抑制半径范围内选取激活（胜利）的column。如果开启了全局抑制，
     * 选举范围为所有columns。
    */
    void inhibitColumns_(UInt activeVector[]);

    /**
     * 返回第columnIndex个column邻居中第k大的overlap值，作为激活时的依据。
    */
    Real kthScore(UInt columnIndex, UInt numActive);

    /**
     * 更新每个column的邻居。在更新抑制半径后需要调用此函数，显式的更新每
     * 个column的邻居。
    */
    void updateNeighbour_();

    /**
     * 调整所有激活column的突触永久值。若突触连接增加永久值，反之减少。
    */
    void adaptSynapses_();

    /**
     * 更新activeDutyCycle_和overlapDutyCycle_。
    */
    void updateDutyCycles_();

    /**
     * 更新所有的BoostFactor。
    */
    void updateBoostFactors_();

    /**
     * 如果overlapDutyCycle(c) < minOverlapDutyCycle(c)那么需要增加该
     * column所有突触的永久值。
    */
    void bumpUpWeakColumns_();

    /**
     * 更新每个column的minOverlapDutyCycle。
    */
    void updateMinDutyCycles_();

    /**
     * 为第colIndex个column所有的突触的永久值增加incFactor。
     * @param colIndex column的索引。
     * @param incFactor 增涨幅度。permenence += incFactor
    */
    void increasePermanences_(UInt colIndex, Real incFactor);

    vector<UInt> getNeighbours_(UInt colIdx);

private:
    UInt numInputs_;
    UInt numColumns_;
    vector<UInt> columnDimensions_;
    vector<UInt> inputDimensions_;
    UInt potentialRadius_;
    Real potentialPct_;
    Real initConnectedPct_;
    bool globalInhibition_;
    int numActiveColumnsPerInhArea_;
    Real localAreaDensity_;
    UInt stimulusThreshold_;
    Real synPermInactiveDec_;
    Real synPermActiveInc_;
    Real synPermConnected_;
    Real minPctOverlapDutyCycles_;
    UInt inhibitionRadius_;
    UInt dutyCyclePeriod_;
    Real boostStrength_;
    /*迭代次数——已有多少输入经过compute*/
    UInt iterationNum_;
    

    vector<vector<UInt>> potentialPools_;
    vector<vector<synapse>> potentialSyn_;
    vector<vector<synapse>> connectedSyn_;

    vector<vector<UInt>> neighbours_;
    vector<UInt> overlap_;
    vector<Real> boostFactors_;
    vector<Real> boostedOverlaps_;
    vector<UInt> activeColumns_;

    vector<Real> activeDutyCycle_;
    vector<Real> overlapDutyCycle_;
    vector<Real> minOverlapDutyCycles_;
    vector<UInt> activeCount_;
    vector<UInt> overlapCount_;
};