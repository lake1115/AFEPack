/**
 * @file   MemoryReclaimer.h
 * @author Robert Lie
 * @date   Sun Apr 29 15:48:52 2007
 * 
 * @brief  
 * 
 * 
 */

#ifndef _MemoryReclaimer_h_
#define _MemoryReclaimer_h_

#include "HGeometry.h"

/*!

  ��HGeometryTree��IrregularMesh�У���һ��HGeometry��HElement������ϸ
  �����Ժ����ǽ������洢���ڴ��У�����Щ������ʱ��ʹ�õ�ʱ������
  ��û�����̽���Щ�ڴ��ͷŵ�������Ҫ�Ŀ�����Ϊ������Щ�����ٴ���Ҫ��ʱ
  �򣬿��Խ�ʡ�����ʱ�䡣����������ض���ʱ�������������Щ��ʱ����
  ���ڴ��ջأ����Ա�����ʵ�֣����͵Ĵ��������£�

  <pre>

    HGeometryTree<DIM,DOW> h_tree;
    ... ...
    IrregularMesh<DIM,DOW> ir_mesh_0;
    IrregularMesh<DIM,DOW> ir_mesh_1;
    ... ...
    MemoryReclaimer<DIM,DOW> mr(h_tree);
    mr.addIrregularMesh(ir_mesh_0);
    mr.addIrregularMesh(ir_mesh_1);
    mr.reclaim();

  </pre>

  �ڵ���reclaim����֮�󣬲�ʹ�õ��ڴ潫�ᱻ���ա�

  ��Ҫ�ر�ע����ǣ������˲���ʱ�����н�����h_tree�ϵġ�����󻹻�ʹ��
  ��IrregularMesh������Ҫʹ��addIrregularMesh �������뵽������յĹ���
  �У�������ܳ���Ǳ�ڵĴ���

 */
template <int DIM, int DOW=DIM>
  class MemoryReclaimer {
 public:
 enum {dim = DIM, dow = DOW};

 typedef HGeometryTree<DIM,DOW> tree_t;
 typedef IrregularMesh<DIM,DOW> ir_mesh_t;
 private:
 tree_t * h_tree;
 std::list<ir_mesh_t *> ir_mesh;

 public:
 MemoryReclaimer() : h_tree(NULL) {}
 MemoryReclaimer(tree_t& _h_tree) : h_tree(&_h_tree) {}
 virtual ~MemoryReclaimer() {}

 public:
 /**
  * ���ý����ڴ���ղ�����HGeometryTree��
  * 
  */
 void setGeometryTree(tree_t& _h_tree) {
   h_tree = &_h_tree;
 }
 /**
  * �ӽ�һ���μӲ�����IrregularMesh��
  * 
  */
 void addIrregularMesh(ir_mesh_t& _ir_mesh) {
   if (&(_ir_mesh.geometryTree()) != h_tree) {
     std::cout << "warning: the irregular mesh added is not based on the geometry tree used."
     << std::endl;
   }
   ir_mesh.push_back(&_ir_mesh);
 }
 void clear() {
   h_tree = NULL;
   ir_mesh.clear();
 }
 /*!

   HGeometryTree�м�����м������໥��������Ϊһ���ǳ����ӵ�����ṹ��
   Ϊ�˻����ڴ棬����ʹ������Ĳ�����ʵ�֣�

   1. ��ÿ��IrregularMesh�в��õ��ڴ���л��ա���������Ƿǳ��򵥵ģ�
   ��Ϊ����ֻ��Ҫ��ÿ��Ҷ�ӽڵ������ͨͨɾ���͹��ˣ�ʵ���ں���
   reclaimIrregularMesh��������С�

   2. ������HGeometryTree�е����м�����ʹ��index=-1����ʶ�����������
   initialTreeLabel�н��С�

   3. �����е�IrregularMesh���������ĵ�Ԫ�õ��ļ����嶼���±�ʶΪ1��
   ��������� labelIrregularMesh ��������

   4. ������HGeometryTree�еļ�������б��������һ��������ı�ʶΪ-1��
   ��ü�����Ӧ�ñ�ɾ����������������һ�α���������ʱ�����Ǿͽ�
   ���ʶ�޸�Ϊ-2��������-1��ʾ�ǵ�һ�����������������ʶΪ1����ô
   ������������ڱ�ĳ��IrregularMeshʹ�ã����ǲ���ʲô�����һ����
   �����ʶΪ-2����ô����������Ѿ������ٵڶ��α��������ˣ���ô����
   �϶�����������������ϵ������������������Ǹ�ָ����Ϊ NULL����
   �����������HGeometryTree��ԭ���Ǹ���״�ṹ�����ݣ��Ǹ���Ҫ��ɾ
   ���Ĳ��֣��Ѿ������Ϊ����״�ṹ��

   5. ������HGeometryTree�еļ�������б��������һ��������ı�ʶΪ-2��
   ���Ǿͽ���ɾ����4��5���������� reclaimTreeMemory ��ʵ�ֵģ�

 */
 void reclaim();

 private:
 void reclaimIrregularMesh(ir_mesh_t&);
 void initialTreeLabel();
 void labelIrregularMesh(ir_mesh_t&);
 void reclaimTreeMemory();

 template <int DIM1> void labelHGeometry(HGeometry<DIM1,DOW>&, int lab);
 template <int DIM1> void labelHGeometryRecursively(HGeometry<DIM1,DOW>& g, int lab);
 template <int DIM1> int relabelHGeometryRecursively(HGeometry<DIM1,DOW>& g);
 template <int DIM1> int reclaimHGeometryRecursively(HGeometry<DIM1,DOW>& g);

 void labelHGeometry(HGeometry<0,DOW>&, int lab);
 void labelHGeometryRecursively(HGeometry<0,DOW>& g, int lab);
 int relabelHGeometryRecursively(HGeometry<0,DOW>& g);
 int reclaimHGeometryRecursively(HGeometry<0,DOW>& g);

 virtual void 
 reclaimHGeometry(void * p_geo, int dim) const {
   switch (dim) {
   case 0: delete ((HGeometry<0,DOW> *)(p_geo)); break;
   case 1: delete ((HGeometry<1,DOW> *)(p_geo)); break;
   case 2: delete ((HGeometry<2,DOW> *)(p_geo)); break;
   case 3: delete ((HGeometry<3,DOW> *)(p_geo)); break;
   }
 }

 protected:
 tree_t * get_tree_ptr() const { return h_tree; }
};

#include "MemoryReclaimer.templates.h"

#endif // _MemoryReclaimer_h_

/**
 * end of file
 * 
 */
