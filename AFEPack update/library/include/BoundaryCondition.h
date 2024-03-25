/**
 * @file   BoundaryCondition.h
 * @author Robert Lie
 * @date   Mon Nov 27 12:17:41 2006
 * 
 * @brief  
 * 
 * 
 */

#ifndef __BoundaryCondition_h__
#define __BoundaryCondition_h__

#include <cstdarg>
#include <vector>
#include <set>
#include <map>

#include <lac/sparse_matrix.h>
#include <lac/vector.h>

#include "Miscellaneous.h"

/**
 * �߽������Ļ��ࡣ�����ӵ��һ���߽��������� _type�����Ա���������
 * �ı߽����������ͣ�Ŀǰ�����������ʵ����Ҫ�Ǿ����������壬������
 * ���Զ�����Ļ���ֻ�е��ϱ߽�������
 *
 * ���⣬����໹������һ�����ϱ�ʶ�������ʾ������Щ���ϱ�ʶ������
 * �ȶ�����ʹ�ñ��߽�����������������ͨ�� add_mark ���������й���
 * ע�⣺�߽������Ĳ��ϱ�ʶֻ�ܼӲ��ܼ��ģ�
 *
 * �����麯�� value �� gradient ���������������ӿڣ�һ������߽纯��
 * �ĺ���ֵ��һ������߽纯�����ݶ�������
 *
 */
class BCondition {
 public:
  /// ���ֻ����ı߽����ͣ���Ϊ����ľ�̬����
  static const int DIRICHLET;
  static const int NEUMANN;
  static const int ROBIN;

 public:
  /*@\{*/
  /// ���캯������������
  BCondition() {}
  BCondition(int type) : _type(type) {}
  BCondition(const BCondition& b) :
    _type(b._type), _bmark(b._bmark) {}
  virtual ~BCondition() {}
  /*@\}*/

 public:
  /// ����������
  BCondition& operator=(const BCondition& b) {
    _type = b._type;
    _bmark = b._bmark;

    return *this;
  }

  /*@\{*/
  /// ��д�߽������
  int type() const {return _type;}
  int& type() {return _type;}
  /*@\}*/

  /*@\{*/
  /// ��д���ϱ�ʶ����
  const std::set<int, std::less<int> >& bound_mark() const {return _bmark;}
  std::set<int, std::less<int> >& bound_mark() {return _bmark;}
  /*@\}*/

  /*@\{*/
  /// �����ϱ�ʶ�����м�����ϱ�ʶ
  void add_mark(int bm) {add_one_mark(bm);}
  void add_mark(const std::vector<int>& bm) {
    u_int n = bm.size();
    for (u_int i = 0;i < n;++ i) {
      add_one_mark(bm[i]);
    }
  }
  template <class V>
    void add_mark(u_int n, const V& bm) {
    for (u_int i = 0;i < n;i ++) {
      add_one_mark(bm[i]);
    }
  }
  void add_mark(u_int n, int bm0, int bm1, ...) {
    add_one_mark(bm0);
    add_one_mark(bm1);

    va_list ap;
    va_start(ap, bm1);
    for (u_int i = 2;i < n;i ++) {
      add_one_mark(va_arg(ap, int));
    }
    va_end(ap);
  }
  /*@\}*/

 private:
  void add_one_mark(int bm) {
    if (bm == 0) {
      std::cerr << "���棺���ϱ�ʶ 0 ��ʾ�����ڲ������ܼӸ��߽�������"
                << std::endl;
      return;
    }
    _bmark.insert(bm);
  }

 public:
  /// �߽�������������ֵ����
  virtual void value(const void * p, void * v) const {}
  /// �߽��������������ݶȺ���
  virtual void gradient(const void * p, void * g) const {}

 private:
  int _type; /// �߽������
  std::set<int, std::less<int> > _bmark; /// ���ϱ�ʶ����
};



/**
 * ͨ��ʹ�ú���ָ����ʵ�� value �� gradient ���������ı߽�������
 *
 */
template <class P, class V, class G = std::vector<V> >
class BCFunction : public BCondition {
  public:
  typedef void (*val_fun_t)(const P&, V&);
  typedef void (*grad_fun_t)(const P&, G&);

  private:
  static void _default_val(const P& p, V& v) {v = V();}
  static void _default_grad(const P& p, G& g) {g = G();}

  val_fun_t _val_fun; /// ��ֵ����ָ��
  grad_fun_t _grad_fun; /// ���ݶȺ���ָ��

  public:
  BCFunction(val_fun_t vf = &_default_val, 
             grad_fun_t gf = &_default_grad) :
  _val_fun(vf), 
  _grad_fun(gf) {}

  BCFunction(int type, 
             val_fun_t vf = &_default_val, 
             grad_fun_t gf = &_default_grad) : 
  BCondition(type), 
  _val_fun(vf), 
  _grad_fun(gf) {}

  virtual ~BCFunction() {}

  public:
  /*@\{*/
  /// ��д��������ָ�������
  const val_fun_t& value_fun_ptr() const {
    return _val_fun;
  }
  val_fun_t& value_fun_ptr() {
    return _val_fun;
  }

  const grad_fun_t& gradient_fun_ptr() const {
    return _grad_fun;
  }
  grad_fun_t& gradient_fun_ptr() {
    return _grad_fun;
  }
  /*@\}*/

  public:
  /*@\{*/
  /// ͨ�����ú���ָ��ı�����������ֵ�����ݶȡ�
  virtual void value(const void * p, void * v) const {
    (*_val_fun)(*(const P *)p, *(V *)v);
  }
  virtual void gradient(const void * p, void * g) const {
    (*_grad_fun)(*(const P *)p, *(G *)g);
  }
  /*@\}*/
};


/**
 * �߽�����������������������һ�ѵı߽�������Ȼ���ھ����ض��Ĳ���
 * ��ʶ�����ɶȴ�ʹ����Ӧ�ı߽�������Ŀǰ������໹�ֲܴڣ�ֻ���Զ�
 * ����һ������µĵ��ϱ�ֵ�����һ�����ָ���ǣ�
 *
 *   - ��Ԫ�ϵĻ����������ǵĲ�ֵ�����������ģ�
 *   - ���������Լ���ֵ���ϵ�ֵȡ 1��
 *
 * ����������£�����Ĵ��������������Ľ����
 * 
 * ʹ�ñ������ж�����⼸����ĵ��ʹ��������£�
 *
 * <pre>
 *
 *    BCFunction<Point<DIM>, double> bc(BCondition::DIRICHLET, u_b);
 *    bc.add_mark(5, 1, 2, 3, 4, 5);
 *    BCAdmin bc_admin;
 *    bc_admin.add(bc);
 *    bc_admin.apply(fem_space, mat, u_h, rhs);
 *
 * </pre>
 *
 * ���ж�����һ�����ϱ�ֵ��ȡֵ��ʽ�ɺ��� u_b ������Ȼ������� 5 ����
 * ���ʶ���ֱ�Ϊ 1, 2, 3, 4, 5�������ֵ���뵽�������У�Ȼ��Ӧ�õ���
 * ��������ϡ�
 *
 */
class BCAdmin {
 public:
  typedef BCondition * bc_ptr_t;

 private:
  std::map<int, bc_ptr_t, std::less<int> > _map;

 public:
  /**
   * Ӧ�õ��ϱ�ֵ������һ������ϵͳ�ϡ�
   * 
   */
  template <
    class SP,    /// ����Ԫ�ռ������
    class MAT,   /// ϡ���������ͣ���������Ϊ SparseMatrix<double>
    class XVEC,  /// �����������ͣ���������Ϊ Vector<double>
    class BVEC   /// �Ҷ����������ͣ���������Ϊ Vector<double>
    >
    void apply(const SP& sp,
               MAT& A, 
               XVEC& u, 
               BVEC& f, 
               bool preserve_symmetry = true) 
    {
      u_int n_dof = sp.n_dof();
      const SparsityPattern& spA = A.get_sparsity_pattern();
      const std::size_t * rowstart = spA.get_rowstart_indices();
      const u_int * colnum = spA.get_column_numbers();
      for (u_int i = 0;i < n_dof;++ i) {
        int bm = sp.dofInfo(i).boundary_mark;
        if (bm == 0) continue; /// 0 ȱʡָ�����ڲ�
        const bc_ptr_t bc = this->find(bm);
        /// �������߽��������ɱ����������߲��ǵ��ϱ߽磬������ȥ��
        if (bc == NULL) continue;
        if (bc->type() != BCondition::DIRICHLET) continue;
        bc->value((const void *)(&sp.dofInfo(i).interp_point), 
                  (void *)(&u(i)));

        f(i) = A.diag_element(i)*u(i);
        for (u_int j = rowstart[i] + 1;j < rowstart[i + 1];++ j) {
          A.global_entry(j) -= A.global_entry(j);
        }
        if (preserve_symmetry) {
          for (u_int j = rowstart[i] + 1;j < rowstart[i + 1];++ j) {
            u_int k = colnum[j];
            const u_int * p = std::find(&colnum[rowstart[k] + 1],
                                        &colnum[rowstart[k + 1]], i);
            if (p != &colnum[rowstart[k+1]]) {
              u_int l = p - &colnum[rowstart[0]];
              f(k) -= A.global_entry(l)*f(i)/A.diag_element(i);
              A.global_entry(l) -= A.global_entry(l);
            }
          }
        }
      }
    }

  /**
   * ������ f ���ɱ�������������ɶ���ͬ��ָ���ϵ�Ԫ�����㡣
   * 
   */
  template <
    class SP,  /// ����Ԫ�ռ������
    class VEC  /// ���������ͣ���������Ϊ Vector<double>
    >
    void clear_entry(const SP& sp,
                     VEC& f)
    {
      u_int n_dof = sp.n_dof();
      for (u_int i = 0;i < n_dof;++ i) {
        int bm = sp.dof_info(i).bound_mark();
        if (bm == 0) continue; /// 0 ȱʡָ�����ڲ�
        const bc_ptr_t bc = find(bm);
        if (bc != NULL) f(i) -= f(i);
      }
    }

  /**
   * ���������м���һ���߽�������һ���߽����������뵽���������Ժ�
   * ������������ add_mark �Ժ��ǲ����Զ������õġ�����д�Ҫ��
   * ��Ҫ���µ��ñ�������
   * 
   */
  template <class BC> 
    void add(BC& b) 
    {
      const std::set<int, std::less<int> >& bm = b.bound_mark();
      std::set<int, std::less<int> >::const_iterator 
        the_bm = bm.begin(), end_bm = bm.end();
      for (;the_bm != end_bm;++ the_bm) {
        _map[*the_bm] = &b;
      }
    }

  /**
   * �ҳ���Ӧ�ڲ��ϱ�ʶ bm �ı߽�������ָ�룬���û���ҵ����򷵻�
   * NULL��
   * 
   */
  bc_ptr_t find(int bm) const 
    {
      std::map<int, bc_ptr_t, std::less<int> >::const_iterator
        the_ptr = _map.find(bm);
      if (the_ptr != _map.end()) {
        return the_ptr->second;
      } else {
        return NULL;
      }
    }
};

#endif // __BoundaryCondition_h__

/**
 * end of file
 * 
 */
