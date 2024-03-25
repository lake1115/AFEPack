/**
 * @file   MemoryReclaimer.templates.h
 * @author Robert Lie
 * @date   Sun Apr 29 15:50:15 2007
 * 
 * @brief  �ڴ���յ�ʵ��
 * 
 * 
 */

#define TEMPLATE template <int DIM, int DOW>
#define THIS MemoryReclaimer<DIM,DOW>

/**
 * ������ g ʹ�� lab ���б�ʶ��
 */
TEMPLATE
template <int DIM1> inline void 
THIS::labelHGeometry(HGeometry<DIM1,DOW>& g, int lab)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    labelHGeometry(*(g.vertex[i]), lab);
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    labelHGeometry(*(g.boundary[i]), lab);
  }
  g.index = lab;
}

/**
 * �Լ����� g ���еݹ��ʶ��������к�����ȱ�ʶ������Ȼ���ʶ�Լ���
 * 
 */
TEMPLATE
template <int DIM1> inline void 
THIS::labelHGeometryRecursively(HGeometry<DIM1,DOW>& g, int lab)
{
  for (int i = 0;i < g.n_boundary;++ i) {
    labelHGeometryRecursively(*(g.boundary[i]), lab);
  }
  if (g.isRefined()) {
    for (int i = 0;i < g.n_child;++ i) {
      labelHGeometryRecursively(*(g.child[i]), lab);
    }
  }

  labelHGeometry(g, lab);
}

/**
 * �Ե㼸���� g ʹ�� lab ���б�ʶ��
 * 
 */
TEMPLATE inline void 
THIS::labelHGeometry(HGeometry<0,DOW>& g, int lab)
{
  g.index = lab;
}

TEMPLATE inline void 
THIS::labelHGeometryRecursively(HGeometry<0,DOW>& g, int lab)
{
  g.index = lab;
}

/**
 * ��IrregularMesh��δʹ�õ��ڴ������ա�����ֻ�Ǽ򵥵ؽ�ÿ��Ҷ�ӽڵ��
 * ���к��ɾ�����ɡ�
 * 
 */
TEMPLATE void 
THIS::reclaimIrregularMesh(typename THIS::ir_mesh_t& m) 
{
  ActiveElementIterator<DIM, DOW> 
    the_ele = m.beginActiveElement(),
    end_ele = m.endActiveElement();
  for (;the_ele != end_ele;++ the_ele) {
    if (the_ele->isRefined()) {
      for (int i = 0;i < the_ele->n_child;++ i) {
        m.deleteTree(the_ele->child[i]);
        the_ele->child[i] = NULL;
      }
    }
  }
}

/**
 * ������HGeometryTree�е����м����嶼ʹ�� -1 ���б�ʶ��
 * 
 */
TEMPLATE
void THIS::initialTreeLabel()
{
  typename HGeometryTree<DIM,DOW>::RootIterator 
    the_ele = h_tree->beginRootElement(),
    end_ele = h_tree->endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    labelHGeometryRecursively(*the_ele, -1);
  }
}

/**
 * ����IrregularMesh���������ĵ�Ԫ�������������ʹ�� 1 ���б�ʶ��
 * 
 */
TEMPLATE
void THIS::labelIrregularMesh(typename THIS::ir_mesh_t& m)
{
  RootFirstElementIterator<DIM, DOW>
    the_ele = m.beginRootFirstElement(),
    end_ele = m.endRootFirstElement();
  for (;the_ele != end_ele;++ the_ele) {
    labelHGeometry(*(the_ele->h_element), 1);
  }

  typename THIS::ir_mesh_t::mesh_t& mesh = m.regularMesh();
#define LABEL_USED_HGEOMETRY(D)                                         \
  if (dim >= D) {                                                       \
    int n_geo = mesh.h_geometry()[D].size();                            \
    for (int i = 0;i < n_geo;++ i) {                                    \
      HGeometry<D,dow> * p_geo = mesh.template h_geometry<D>(i);        \
      labelHGeometry(*p_geo, 1);                                        \
    }                                                                   \
  }

  LABEL_USED_HGEOMETRY(0);
  LABEL_USED_HGEOMETRY(1);
  LABEL_USED_HGEOMETRY(2);
  LABEL_USED_HGEOMETRY(3);

#undef LABEL_USED_HGEOMETRY
}

/**
 * �Լ����� g ����������ͺ�������ر�ʶ�������ԭ��ʶΪ 1��������ʲ
 * ô������������ 1�������ԭ��ʶΪ -1��˵������������ǵ�һ�α���������
 * ���ǽ����ʶ�޸�Ϊ -2�������� -1�������ԭ��ʶΪ -2��˵�����������
 * �Ѿ����ǵ�һ�α��������ˣ����ǽ�������������Ȼ�󷵻� -2��
 * 
 */
TEMPLATE
template <int DIM1> inline int 
THIS::relabelHGeometryRecursively(HGeometry<DIM1,DOW>& g)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    if (g.vertex[i] == NULL) continue;
    if (relabelHGeometryRecursively(*(g.vertex[i])) == -2) {
      g.vertex[i] = NULL;
    }
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    if (g.boundary[i] == NULL) continue;
    if (relabelHGeometryRecursively(*(g.boundary[i])) == -2) {
      g.boundary[i] = NULL;
    }
  }
  if (g.isRefined()) {
    for (int i = 0;i < g.n_child;++ i) {
      if (g.child[i] == NULL) continue;
      if (relabelHGeometryRecursively(*(g.child[i])) == -2) {
        g.child[i] = NULL;
      }
    }
  }
  Assert ((g.index == -1) || 
          (g.index ==  1) ||
          (g.index == -2), ExcInternalError());
  if (g.index == -1) {
    g.index = -2;
    return -1;
  }
  else return g.index;
}


/**
 * �Ե㼸���� g �������±�ʶ��
 * 
 */
TEMPLATE inline int 
THIS::relabelHGeometryRecursively(HGeometry<0,DOW>& g)
{
  Assert ((g.index == -1) || 
          (g.index ==  1) ||
          (g.index == -2), ExcInternalError());
  if (g.index == -1) {
    g.index = -2;
    return -1;
  }
  else return g.index;
}

/**
 * �Լ����� g ��������Ѿ��������л��ա�
 * 
 */
TEMPLATE
template <int DIM1> inline int 
THIS::reclaimHGeometryRecursively(HGeometry<DIM1,DOW>& g)
{
  for (int i = 0;i < g.n_vertex;++ i) {
    if (g.vertex[i] != NULL) {
      if (reclaimHGeometryRecursively(*(g.vertex[i])) == -1) {
        g.vertex[i] = NULL;
      }
    }
  }
  for (int i = 0;i < g.n_boundary;++ i) {
    if (g.boundary[i] != NULL) {
      if (reclaimHGeometryRecursively(*(g.boundary[i])) == -1) {
        g.boundary[i] = NULL;
      }
    }
  }
  for (int i = 0;i < g.n_child;++ i) {
    if (g.child[i] == NULL) continue;
    if (reclaimHGeometryRecursively(*(g.child[i])) == -1) {
      g.child[i] = NULL;
    }
  }
  Assert ((g.index == 1 || g.index == -2), ExcInternalError());
  if (g.index == -2) {
    this->reclaimHGeometry(&g, DIM1);
    return -1;
  }
  else {
    return 1;
  }
}


/**
 * �Ե㼸���� g ���л��ա�
 * 
 */
TEMPLATE inline int 
THIS::reclaimHGeometryRecursively(HGeometry<0,DOW>& g)
{
  Assert ((g.index == 1 || g.index == -2), ExcInternalError());
  if (g.index == -2) {
    this->reclaimHGeometry(&g, 0);
    return -1;
  }
  else {
    return 1;
  }
}


/**
 * ���� HGeometryTree �в��õ��ڴ档
 * 
 */
TEMPLATE
void THIS::reclaimTreeMemory()
{
  typename HGeometryTree<DIM,DOW>::RootIterator
    the_ele = h_tree->beginRootElement(),
    end_ele = h_tree->endRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    relabelHGeometryRecursively(*the_ele);
  }

  the_ele = h_tree->beginRootElement();
  for (;the_ele != end_ele;++ the_ele) {
    reclaimHGeometryRecursively(*the_ele);
  }
}

TEMPLATE
void THIS::reclaim()
{
  /// �ȶԸ���IrregularMesh�Ĳ�ʹ���ڴ���л���
  typename std::list<typename THIS::ir_mesh_t *>::iterator 
    the_ir_mesh = ir_mesh.begin(),
    end_ir_mesh = ir_mesh.end();
  for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
    reclaimIrregularMesh(**the_ir_mesh);
  }
  
  /// �����������־Ϊ -1
  initialTreeLabel();

  /// ��ʹ���е������ʶΪ 1
  the_ir_mesh = ir_mesh.begin();
  for (;the_ir_mesh != end_ir_mesh;++ the_ir_mesh) {
    labelIrregularMesh(**the_ir_mesh);
  }

  /// ���ڴ���л���
  reclaimTreeMemory();
}

#undef THIS
#undef TEMPLATE

/**
 * end of file
 * 
 */
