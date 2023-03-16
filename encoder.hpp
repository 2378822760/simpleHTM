/** @file
 * ���ļ�����ϵͳ��ʹ�õı�����������ͬ���͵��������ݱ����HTM���õ�ϡ��ֲ�������SDR����
 * ����������ȷ����Щ�������Ӧ���� 1����ЩӦ���� 0���Ӷ��������ݵ���Ҫ����������
 * ���Ƶ�������ϢӦ�ò����˴˸߶��ص��� SDR��
 */
#include "type.hpp"

#include <vector>

 /** ����������
  * �ڱ�������ʱ����һЩ��Ҫ�ķ�����Ҫ���ǣ�
  * 1. �������Ƶ�����Ӧ�ò�����Ծ���ر˴��ص��� SDR
  * 2. ��ͬ����������Ӧ�����ǲ�����ͬ��������
  * 3. �����е�������Ϣ�������ϢӦ������ͬ��ά�ȣ�����������
  * 4. �����е�������Ϣ�������ϢӦ�������Ƶ�ϡ��ȣ��������㹻���� 1 ���ش����������Ӳ���
  */
class ScalarEncoder {
public:
    /**
     * ����һ������������
     *
     * @param w ���뵥����ֵ��Ҫ�ı����� -- ����źŵġ�width����һ�㲻����20~25��
     * @param minValue ����������
     * @param maxValue ����������
     * @param bucketNum Ͱ������
     * @param clipInput �����볬�����䣬��Ϊtrue�ض���ֵ������������
     */
    ScalarEncoder(int w, double minValue, double maxValue, int bucketNum, bool clipInput);


    /**
     * ��ȡ������ݳ���
     */
    UInt getOutputWidth() const {
        return n_;
    }
    /**
     * ��������ֵ�ı�����д��output
     * output��Ҫ��ǰ������㹻�Ŀռ䣨w��
     */
    void encodeIntoArray(double input, UInt output[]);

private:
    int w_;
    double minValue_;
    double maxValue_;
    int bucketNum_;
    UInt n_;
    bool clipInput_;
}; // end class ScalarEncoder