/**
 * @file   SparseMatrixTool.h
 * @author Robert Lie
 * @date   Sun Nov  4 09:06:11 2007
 * 
 * @brief  
 * 
 * 
 */

#include <base/exceptions.h>
#include <lac/sparse_matrix.h>

#include "Miscellaneous.h"

/**
 * ����һЩ��ϡ������ϡ�����ģ����в����Ĺ��ߺ���
 * 
 */
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
			   SparsityPattern& sp); 

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
			   SparsityPattern& sp);

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
   * ��ģ�塣Ҫ�� \f$ A \f$ �� \f$ B \f$ ������ͬ������
   * �Լ� \f$ A \f$ ��һ������
   * 
   * @param sp0 ���� \f$ A \f$ ��ģ��
   * @param sp1 ���� \f$ B \f$ ��ģ��
   * @param sp �õ�����ģ��
   */
  void gammaCatSparsityPattern(const SparsityPattern& sp0,
			       const SparsityPattern& sp1,
			       SparsityPattern& sp);

  /**
   * ������ϡ�����ģ�尴�նԽǷ�ʽ���ӳ�Ϊһ���µ�ģ�塣
   * 
   * @param sp0 ���ϽǾ����ģ��
   * @param sp1 ���½Ǿ����ģ��
   * @param sp �õ�����ģ��
   */
  void dCatSparsityPattern(const SparsityPattern& sp0,
			   const SparsityPattern& sp1,
			   SparsityPattern& sp);

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
			      SparsityPattern& sp);

  /**
   * ������ϡ�������ˮƽ��ʽ���ӳ�Ϊһ���µ�ϡ�����
   * 
   * @param m0 ��ߵľ���
   * @param m1 �ұߵľ���
   * @param sp �õ����¾����ģ��
   * @param m �õ����¾���
   * @param is_pattern_ok ���Ϊ�棬���ģ����й��죻����
   *                      �ٶ�ģ���Ѿ���������ȷ�Ĺ��죻ȱ
   *                      ʡֵΪ�档
   */
  template <typename number>
    void hCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok = true);

  /**
   * ������ϡ���������ֱ��ʽ���ӳ�Ϊһ���µ�ϡ�����
   * 
   * @param m0 �ϱߵľ���
   * @param m1 �±ߵľ���
   * @param sp �õ����¾����ģ��
   * @param m �õ����¾���
   * @param is_pattern_ok ���Ϊ�棬���ģ����й��죻����
   *                      �ٶ�ģ���Ѿ���������ȷ�Ĺ��죻ȱ
   *                      ʡֵΪ�档
   */
  template <typename number>
    void vCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok = true);

  /**
   * ������ϡ����󰴰�������ķ�ʽ���ӳ�Ϊһ���µ�ϡ��
   * ����Ҳ���Ƕ��ڸ�����ϡ����� \f$ A \f$ �� \f$ B
   * \f$�����ǵõ�����
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A & B \\
   *      B^T & 0
   *    \end{array}\right)
   * \f]
   *
   * @param m0 ���� \f$ A \f$
   * @param m1 ���� \f$ B \f$
   * @param sp �õ����¾����ģ��
   * @param m �õ����¾���
   * @param is_pattern_ok ���Ϊ�棬���ģ����й��죻����
   *                      �ٶ�ģ���Ѿ���������ȷ�Ĺ��죻ȱ
   *                      ʡֵΪ�档
   */
  template <typename number>
    void gammaCatSparseMatrix(const SparseMatrix<number>& m0,
			      const SparseMatrix<number>& m1,
			      SparsityPattern& sp,
			      SparseMatrix<number>& m,
			      bool is_pattern_ok = true);

  /**
   * ������ϡ������նԽǷ�ʽ���ӳ�Ϊһ���µ�ϡ�����
   * 
   * @param m0 ���ϽǾ���
   * @param m1 ���½Ǿ���
   * @param sp �õ����¾����ģ��
   * @param m �õ����¾���
   * @param is_pattern_ok ���Ϊ�棬���ģ����й��죻����
   *                      �ٶ�ģ���Ѿ���������ȷ�Ĺ��죻ȱ
   *                      ʡֵΪ�档
   */
  template <typename number>
    void dCatSparseMatrix(const SparseMatrix<number>& m0,
			  const SparseMatrix<number>& m1,
			  SparsityPattern& sp,
			  SparseMatrix<number>& m,
			  bool is_pattern_ok = true);

  /** 
   * ���ĸ�ϡ�������ȫ����Ϊһ���µ�ϡ�����
   *
   * \f[
   *    \left(\begin{array}{cc}
   *      A_{00} & A_{01} \\
   *      A_{10} & A_{11}
   *    \end{array}\right)
   * \f]
   * 
   * @param m00 ���ϽǾ���
   * @param m01 ���ϽǾ���
   * @param m10 ���½Ǿ���
   * @param m11 ���½Ǿ���
   * @param sp �õ����¾����ģ��
   * @param m �õ����¾���
   * @param is_pattern_ok ���Ϊ�棬���ģ����й��죻����
   *                      �ٶ�ģ���Ѿ���������ȷ�Ĺ��죻ȱ
   *                      ʡֵΪ�档
   */  
  template <typename number>
    void fullCatSparseMatrix(const SparseMatrix<number>& m00,
			     const SparseMatrix<number>& m01,
			     const SparseMatrix<number>& m10,
			     const SparseMatrix<number>& m11,
			     SparsityPattern& sp,
			     SparseMatrix<number>& m,
			     bool is_pattern_ok = true);

  DeclException0(ExcNotCompressed);
  DeclException2(ExcDimensionDontMatch, int, int,
		 << "The dimensions " << arg1 << " and " << arg2
		 << " do not match properly.");
};

/**
 * end of file
 * 
 */
