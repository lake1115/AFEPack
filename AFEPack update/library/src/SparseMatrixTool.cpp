/**
 * SparseMatrixTool.cpp
 * 
 */

#include "SparseMatrixTool.templates.h"

namespace SparseMatrixTool {

  /**
   * ������ϡ�����ģ�尴ˮƽ�������ӳ�Ϊһ���µ�ģ��
   * Ҫ�����������������ͬ��������
   * 
   * @param sp0 ��ߵ�ģ��
   * @param sp1 �ұߵ�ģ��
   * @param sp �õ�����ģ��
   */
  void hCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows();
    u_int n0 = sp0.n_cols();

    /**
     * ����������ģ���Ƿ��Ѿ�ѹ�����ˡ�
     * 
     */
    Assert(sp0.is_compressed(), SparsityPattern::ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * ���������ģ���Ƿ������ͬ��������
     * 
     */
    Assert (m0 == sp1.n_rows(), ExcDimensionDontMatch(m0, sp1.n_rows()));
  
    int i, j;
    u_int n1 = sp1.n_cols();
    std::vector<u_int> row_length(m0, 0);
    for (i = 0;i < m0;i ++)
      row_length[i] = sp0.row_length(i) + sp1.row_length(i);
    sp.reinit(m0, n0 + n1, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++) {
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(i, n0 + p_column1[j]);
    }
    sp.compress();
  };

  /**
   * ������ϡ�����ģ�尴��ֱ�������ӳ�Ϊһ���µ�ģ��
   * Ҫ�����������������ͬ��������
   *
   * @param sp0 �ϱߵ�ģ��
   * @param sp1 �±ߵ�ģ��
   * @param sp �õ�����ģ��
   */
  void vCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows();
    u_int n0 = sp0.n_cols();

    /**
     * ����������ģ���Ƿ��Ѿ�ѹ�����ˡ�
     * 
     */
    Assert(sp0.is_compressed(), ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * ���������ģ���Ƿ������ͬ��������
     * 
     */
    Assert (n0 == sp1.n_cols(), ExcDimensionDontMatch(m0, sp1.n_cols()));
  
    u_int i, j, m1 = sp1.n_rows();
    std::vector<u_int> row_length(m0 + m1, 0);
    for (i = 0;i < m0;i ++)
      row_length[i] = sp0.row_length(i);
    for (i = 0;i < m1;i ++)
      row_length[i + m0] = sp1.row_length(i);
    sp.reinit(m0 + m1, n0, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++)
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
    for (i = 0;i < m1;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(m0 + i, p_column1[j]);
    sp.compress();
  };

  /**
   * ������ϡ�����ģ�尴��������ķ�ʽ���ӳ�Ϊһ���µ�
   * ģ�壬Ҳ���Ƕ��ڸ�����ϡ����� \f$ A \f$ �� \f$ B
   * \f$�����ǵõ�����
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A & B \\
   *      B^T & 0
   *    \end{array}\right)
   * \f]
   *
   * ��ģ�塣Ҫ�� \f$ A \f$ �� \f$ B \f$ ������ͬ��������
   * 
   * @param sp0 ���� \f$ A \f$ ��ģ��
   * @param sp1 ���� \f$ B \f$ ��ģ��
   * @param sp �õ�����ģ��
   */
  void gammaCatSparsityPattern(const SparsityPattern& sp0,
			       const SparsityPattern& sp1,
			       SparsityPattern& sp)
  {
    u_int m = sp0.n_rows();

    /**
     * ����������ģ���Ƿ��Ѿ�ѹ�����ˡ�
     * 
     */
    Assert(sp0.is_compressed(), ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    /**
     * ���������ģ���Ƿ������ͬ�������Լ���һ��ģ���Ƿ�
     * ��һ�������ģ�塣
     * 
     */
    Assert (m == sp0.n_cols(), ExcDimensionDontMatch(m, sp0.n_cols()));
    Assert (m == sp1.n_rows(), ExcDimensionDontMatch(m, sp1.n_rows()));

    int i, j;
    u_int n1 = sp1.n_cols();
    std::vector<u_int> row_length(m + n1, 1);
    for (i = 0;i < m;i ++)
      row_length[i] = sp0.row_length(i) + sp1.row_length(i);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	row_length[m + p_column1[j]] += 1;
    sp.reinit(m + n1, m + n1, row_length);
    for (i = 0;i < m;i ++) {
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++) {
	sp.add(i, m + p_column1[j]);
	sp.add(m + p_column1[j], i);
      }
    }
    sp.compress();
  };

  /** 
   * ������ϡ�����ģ�尴�ԽǷ������ӳ�Ϊһ���µ�ģ��
   * 
   * @param sp0 ���Ͻǵ�ģ��
   * @param sp1 ���½ǵ�ģ��
   * @param sp �õ�����ģ��
   */
  void dCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp)
  {
    u_int m0 = sp0.n_rows(), m1 = sp1.n_rows();
    u_int n0 = sp0.n_cols(), n1 = sp1.n_cols();

    /**
     * ����������ģ���Ƿ��Ѿ�ѹ�����ˡ�
     * 
     */
    Assert(sp0.is_compressed(), SparsityPattern::ExcNotCompressed());
    Assert(sp1.is_compressed(), ExcNotCompressed());

    int i = 0, j;
    if (m0 + m1 == n0 + n1) i = 1;
    std::vector<u_int> row_length(m0 + m1, i);
    for (i = 0;i < m0;i ++)
      row_length[i] += sp0.row_length(i); 
    for (i = 0;i < m1;i ++)
      row_length[i + m0] += sp1.row_length(i);
    sp.reinit(m0 + m1, n0 + n1, row_length);
    const std::size_t * p_rowstart0 = sp0.get_rowstart_indices();
    const u_int * p_column0 = sp0.get_column_numbers();
    const std::size_t * p_rowstart1 = sp1.get_rowstart_indices();
    const u_int * p_column1 = sp1.get_column_numbers();
    for (i = 0;i < m0;i ++)
      for (j = p_rowstart0[i];j < p_rowstart0[i + 1];j ++)
	sp.add(i, p_column0[j]);
    for (i = 0;i < m1;i ++)
      for (j = p_rowstart1[i];j < p_rowstart1[i + 1];j ++)
	sp.add(m0 + i, n0 + p_column1[j]);
    sp.compress();
  };


  /** 
   * ���ĸ�ϡ�����ģ����ȫ����Ϊһ���µ�ģ��
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A_{00} & A_{01} \\
   *      A_{01} & A_{11}
   *    \end{array}\right)
   * \f]
   * 
   * @param sp00 ���ϽǾ����ģ��
   * @param sp01 ���ϽǾ����ģ��
   * @param sp10 ���½Ǿ����ģ��
   * @param sp11 ���½Ǿ����ģ��
   * @param sp �õ�����ģ��
   */  
  void fullCatSparsityPattern(const SparsityPattern& sp00,
			      const SparsityPattern& sp01,
			      const SparsityPattern& sp10,
			      const SparsityPattern& sp11,
			      SparsityPattern& sp)
  {
    SparsityPattern sp0, sp1;
    hCatSparsityPattern(sp00, sp01, sp0);
    hCatSparsityPattern(sp10, sp11, sp1);
    vCatSparsityPattern(sp0, sp1, sp);
    sp.compress();
  };

#define number double
  template void hCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);
  template void vCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);
  template void gammaCatSparseMatrix<number>(const SparseMatrix<number>&,
					     const SparseMatrix<number>&,
					     SparsityPattern&,
					     SparseMatrix<number>&,
					     bool);
  template void dCatSparseMatrix<number>(const SparseMatrix<number>&,
					 const SparseMatrix<number>&,
					 SparsityPattern&,
					 SparseMatrix<number>&,
					 bool);

  template void fullCatSparseMatrix(const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    const SparseMatrix<number>&,
				    SparsityPattern&,
				    SparseMatrix<number>&,
				    bool is_pattern_ok);
#undef number
};

/**
 * end of file
 * 
 */
