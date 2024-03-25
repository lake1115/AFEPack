/**
 * @file   ThreadBase.h
 * @author Ruo Li
 * @date   Sun Feb 13 20:08:19 2005
 * 
 * @brief  for convenient using of POSIX thread in C++
 * 
 * 
 */

/*!
\class Thread?
\brief ������Ŀ����Ϊ���ܹ��ȽϷ������AFEPack��ʹ��POSIX�߳��ṩһ���ӿڡ�

������ʵ�ֱȽϸ��ӣ�����Ҫ�ǽ����deal.II�й��ڶ��߳�֧�ֵķ����������ǳ�
���˽����еļ���������ȥ�ο�������ĵ���ѧϰһ�¿������C++�����ˮƽ������
�пյ�ʱ��ȥ��һ�ۡ����������д�ñȽ��Ѷ������鲻Ҫ���ҵ�������룬�˷�ʱ�䡣
���������ˣ�˵�����ˮƽ���Ҹߣ�;-)�����⣬�Ҳ��ṩ�κι����̰߳�ȫ�Է���
�ı�֤��ȫ�������Լ��Ĵ���ȥʵ���ˡ�

deal.II�еĶ��߳�������ACE����������
ֱ���Լ�д��һ���Ƚϼ򵥵��̹߳�������������Ҫ˵��һ��ʹ�õķ���������ڴ��
���û���˵���Ѿ��㹻�ˡ�������̵߳�֧�֣����԰����㷽��Ĳ���һ���̣߳�Эͬ
���һ���Ƚϴ�������̵߳���ڳ�����������Ա������Ҳ���Բ������Ա������
��������������ڳ��������Ա������������÷���

\verbatim
���� �����̵߳���ں����������i��Chebyshev����ʽ�ڵ�0.64����ֵ
void fun(int i, const double& d)
{
  d = cos(i*acos(0.64));
};

int main(int argc, char * argv[])
{
  double val[5];
  ThreadManager th_man;        ���� �����̹߳���������
  for (int i = 0;i < 5;i ++) { ���� ���Ǵ���5���߳�������
    ���� �����̣߳���ʵ�����е����ܶ���������һ����
    th_man.start(encapsulate(&fun).collectArgs(i, val[i])); 
  }
  th_man.join(encapsulate(&fun)); ���� �ȴ��߳̽���
  ���� ������һ��ʵ�ֵ÷ǳ��ѿ���������ʵ����û�а취��
  ���� �������ʲô���У���������θ�����
  for (int i = 0;i < 5;i ++) { ���� ���������
    std::cout << "Chebyshev(" << i 
	      << ", 0.64) = " << val[i]
	      << std::endl;
  }
  return 0;
};
\endverbatim

��������Ա��������ں���������Ȼ��һ���������һ��Ҫ�Ѷ���
��ȥ���������������ӣ�

\verbatim
class Chebyshev {
 private:
  double x;
 public: ���� �̵߳���ں��������ǹ��е�
  Chebyshev(const double& _x) : x(_x) {};
  void fun(int i, const double& d) {
    d = cos(i*acos(x));
  };
};

int main(int argc, char * argv[])
{
  double val[5];
  Chebyshev ch(0.64);          ���� ������Ķ���
  ThreadManager th_man;        ���� �����̹߳���������
  for (int i = 0;i < 5;i ++) { ���� ���Ǵ���5���߳�������
    ���� �����̣߳���������Ҫ�Ѷ��󴫽�ȥ
    th_man.start(encapsulate(&Chebyshev::fun).collectArgs(&ch, i, val[i])); 
  }
  th_man.join(encapsulate(&Chebyshev::fun)); ���� �ȴ��߳̽���
  for (int i = 0;i < 5;i ++) { ���� ���������
    std::cout << "Chebyshev(" << i 
	      << ", 0.64) = " << val[i]
	      << std::endl;
  }
  return 0;
};
\endverbatim

������ʵ�֣�Ψһ��Ŀ�����Ϊ��ʹ���������㡣

����ڽ�������Ԫ�����ʱ����ʵ�ֵ���һ��
\verbatim
setThread(n); ���� n>1
\endverbatim
��ô�ڽ�������Ԫ�ռ�Ĺ����ʱ�����һ���ִ����ʹ�ö��߳̽��У��Ժ�
�һὫ����Ĳ��ֽ��ж��̻߳���������Ҫʱ�䣬�����ĵȴ�Ŷ��

����ʹ��AFEPackд������Ԫ�ĳ��򣬳���������ĺ���
\verbatim
void fun()
{
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    ... ... �������ǻ��������Ԫ�Ͻ���һ���Ƚϸ��ӵ���ֵ���ֵĲ���
  }
};

����������ĵط������ǵ��������������
fun();
\endverbatim
����������һ�γ������ǾͿ���ʹ������ķ�ʽ�����ж��̻߳��ˣ�
\verbatim
void fun(int rank, int n_thread)
{
  int n_element = fem_space.n_elemen();
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
  int stride = n_element/n_thread;
  the_ele += stride*rank;
  if (rank + 1 < n_thread)
    end_ele = the_ele + stride;
  for (;the_ele != end_ele;++ the_ele) { ��������ֻ����һ���ֵ�Ԫ����
    ... ... ���������������Ԫ�Ͻ���һ���Ƚϸ��ӵ���ֵ���ֵĲ���
  }
};

������Ȼ��ʱ���������ĵط���������������΢����һ����
int n_thread = getThread();
ThreadManager th_man;
for (int rank = 1;rank < n_thread;rank ++) { ��������һ���߳�
  th_man.start(encapsulate(&fun).collectArgs(rank, n_thread));
}
fun(0, n_thread); �����Լ�Ҳ��Ҫ����
th_man.join(encapsulate(&fun)); �������Ź������
\endverbatim
�������Լ�����������ʹ�÷�ʽ�����ǱȽϷ���ģ�;-)

*/
template <typename Return
#if ARG>=1 
,typename Arg1 
#endif
#if ARG>=2 
,typename Arg2 
#endif
#if ARG>=3 
,typename Arg3 
#endif
#if ARG>=4 
,typename Arg4 
#endif
#if ARG>=5 
,typename Arg5 
#endif
#if ARG>=6 
,typename Arg6 
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
,typename Class=NullThread>
struct CLASSNAME 
{
  public:

typedef Return (*FunPtr)(
#if ARG>=1 
			 Arg1 
#endif
#if ARG>=2 
			 ,Arg2 
#endif
#if ARG>=3 
			 ,Arg3 
#endif
#if ARG>=4 
			 ,Arg4 
#endif
#if ARG>=5 
			 ,Arg5 
#endif
#if ARG>=6 
			 ,Arg6 
#endif
#if ARG>=7 
			 ,Arg7
#endif
#if ARG>=8 
			 ,Arg8 
#endif
#if ARG>=9 
			 ,Arg9 
#endif
#if ARG>=10 
			 ,Arg10 
#endif
			 );

  struct FunData 
  {
    public:
    FunPtr     fun_ptr;
#if ARG>=1 
    Arg1 arg1;
#endif
#if ARG>=2 
    Arg2 arg2;
#endif
#if ARG>=3 
    Arg3 arg3;
#endif
#if ARG>=4 
    Arg4 arg4;
#endif
#if ARG>=5 
    Arg5 arg5;
#endif
#if ARG>=6 
    Arg6 arg6;
#endif
#if ARG>=7 
    Arg7 arg7;
#endif
#if ARG>=8 
    Arg8 arg8;
#endif
#if ARG>=9 
    Arg9 arg9;
#endif
#if ARG>=10 
    Arg10 arg10;
#endif
    FunData(FunPtr _fun_ptr
#if ARG>=1 
	    ,Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) :
      fun_ptr(_fun_ptr)
#if ARG>=1 
	 ,arg1(_arg1)
#endif
#if ARG>=2 
	 ,arg2(_arg2)
#endif
#if ARG>=3 
	 ,arg3(_arg3)
#endif
#if ARG>=4 
	 ,arg4(_arg4)
#endif
#if ARG>=5 
	 ,arg5(_arg5)
#endif
#if ARG>=6 
	 ,arg6(_arg6)
#endif
#if ARG>=7 
	 ,arg7(_arg7)
#endif
#if ARG>=8 
	 ,arg8(_arg8)
#endif
#if ARG>=9 
	 ,arg9(_arg9)
#endif
#if ARG>=10 
	 ,arg10(_arg10)
#endif
    {};
    FunData(const FunData& _fd) :
      fun_ptr(_fd.fun_ptr)
#if ARG>=1 
	 ,arg1(_fd.arg1)
#endif
#if ARG>=2 
	 ,arg2(_fd.arg2)
#endif
#if ARG>=3 
	 ,arg3(_fd.arg3)
#endif
#if ARG>=4 
	 ,arg4(_fd.arg4)
#endif
#if ARG>=5 
	 ,arg5(_fd.arg5)
#endif
#if ARG>=6 
	 ,arg6(_fd.arg6)
#endif
#if ARG>=7 
	 ,arg7(_fd.arg7)
#endif
#if ARG>=8 
	 ,arg8(_fd.arg8)
#endif
#if ARG>=9 
	 ,arg9(_fd.arg9)
#endif
#if ARG>=10 
	 ,arg10(_fd.arg10)
#endif
    {};
    FunData& operator=(const FunData& _fd) 
{
fun_ptr = _fd.fun_ptr;
#if ARG>=1 
 arg1 = _fd.arg1;
#endif
#if ARG>=2 
 arg2 = _fd.arg2;
#endif
#if ARG>=3 
 arg3 = _fd.arg3;
#endif
#if ARG>=4 
 arg4 = _fd.arg4;
#endif
#if ARG>=5 
 arg5 = _fd.arg5;
#endif
#if ARG>=6 
 arg6 = _fd.arg6;
#endif
#if ARG>=7 
 arg7 = _fd.arg7;
#endif
#if ARG>=8 
 arg8 = _fd.arg8;
#endif
#if ARG>=9 
 arg9 = _fd.arg9;
#endif
#if ARG>=10 
 arg10 = _fd.arg10;
#endif
 return *this;
};
    static void * thread_entry(void *arg_ptr) {
      FunData * mfd = reinterpret_cast<FunData *>(arg_ptr);
      (*(mfd->fun_ptr))(
#if ARG>=1 
			    mfd->arg1
#endif
#if ARG>=2 
			    , mfd->arg2
#endif
#if ARG>=3 
			    , mfd->arg3
#endif
#if ARG>=4 
			    , mfd->arg4
#endif
#if ARG>=5 
			    , mfd->arg5
#endif
#if ARG>=6 
			    , mfd->arg6
#endif
#if ARG>=7 
			    , mfd->arg7
#endif
#if ARG>=8 
			    , mfd->arg8
#endif
#if ARG>=9 
			    , mfd->arg9
#endif
#if ARG>=10 
			    , mfd->arg10
#endif
			  );
    return NULL;
  };

  };

    class ArgCollector
  {
  public:
    typedef FunData FunDataType;
    FunPtr fun_ptr;

    ArgCollector (FunPtr _fun_ptr_) {
      fun_ptr = _fun_ptr_;
    };
        
    FunData collectArgs (
#if ARG>=1 
	    Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) {
      return FunData(fun_ptr
#if ARG>=1 
		     , _arg1
#endif
#if ARG>=2 
		     , _arg2
#endif
#if ARG>=3 
		     , _arg3
#endif
#if ARG>=4 
		     , _arg4
#endif
#if ARG>=5 
		     , _arg5
#endif
#if ARG>=6 
		     , _arg6
#endif
#if ARG>=7 
		     , _arg7
#endif
#if ARG>=8 
		     , _arg8
#endif
#if ARG>=9 
		     , _arg9
#endif
#if ARG>=10 
		     , _arg10
#endif
		      );
    };

  };

typedef Return (Class::*MemFunPtr)(
#if ARG>=1 
				   Arg1 
#endif
#if ARG>=2 
				   ,Arg2 
#endif
#if ARG>=3 
				   ,Arg3 
#endif
#if ARG>=4 
				   ,Arg4 
#endif
#if ARG>=5 
				   ,Arg5 
#endif
#if ARG>=6 
				   ,Arg6 
#endif
#if ARG>=7 
				   ,Arg7
#endif
#if ARG>=8 
				   ,Arg8 
#endif
#if ARG>=9 
				   ,Arg9 
#endif
#if ARG>=10 
				   ,Arg10 
#endif
				   );

  struct MemFunData 
  {
    public:
    MemFunPtr     mem_fun_ptr;
    Class *       object;
#if ARG>=1 
    Arg1 arg1;
#endif
#if ARG>=2 
    Arg2 arg2;
#endif
#if ARG>=3 
    Arg3 arg3;
#endif
#if ARG>=4 
    Arg4 arg4;
#endif
#if ARG>=5 
    Arg5 arg5;
#endif
#if ARG>=6 
    Arg6 arg6;
#endif
#if ARG>=7 
    Arg7 arg7;
#endif
#if ARG>=8 
    Arg8 arg8;
#endif
#if ARG>=9 
    Arg9 arg9;
#endif
#if ARG>=10 
    Arg10 arg10;
#endif
    MemFunData(MemFunPtr _mem_fun_ptr
	       ,Class * _object
#if ARG>=1 
	       ,Arg1 _arg1
#endif
#if ARG>=2 
	       ,Arg2 _arg2
#endif
#if ARG>=3 
	       ,Arg3 _arg3
#endif
#if ARG>=4 
	       ,Arg4 _arg4
#endif
#if ARG>=5 
	       ,Arg5 _arg5
#endif
#if ARG>=6 
	       ,Arg6 _arg6
#endif
#if ARG>=7 
	       ,Arg7 _arg7
#endif
#if ARG>=8 
	       ,Arg8 _arg8
#endif
#if ARG>=9 
	       ,Arg9 _arg9
#endif
#if ARG>=10 
	       ,Arg10 _arg10
#endif
	       ) :
      mem_fun_ptr(_mem_fun_ptr)
	 ,object(_object)
#if ARG>=1 
	 ,arg1(_arg1)
#endif
#if ARG>=2 
	 ,arg2(_arg2)
#endif
#if ARG>=3 
	 ,arg3(_arg3)
#endif
#if ARG>=4 
	 ,arg4(_arg4)
#endif
#if ARG>=5 
	 ,arg5(_arg5)
#endif
#if ARG>=6 
	 ,arg6(_arg6)
#endif
#if ARG>=7 
	 ,arg7(_arg7)
#endif
#if ARG>=8 
	 ,arg8(_arg8)
#endif
#if ARG>=9 
	 ,arg9(_arg9)
#endif
#if ARG>=10 
	 ,arg10(_arg10)
#endif
    {};
    MemFunData(const MemFunData& _mfd) :
      mem_fun_ptr(_mfd.mem_fun_ptr)
	 ,object(_mfd.object)
#if ARG>=1 
	 ,arg1(_mfd.arg1)
#endif
#if ARG>=2 
	 ,arg2(_mfd.arg2)
#endif
#if ARG>=3 
	 ,arg3(_mfd.arg3)
#endif
#if ARG>=4 
	 ,arg4(_mfd.arg4)
#endif
#if ARG>=5 
	 ,arg5(_mfd.arg5)
#endif
#if ARG>=6 
	 ,arg6(_mfd.arg6)
#endif
#if ARG>=7 
	 ,arg7(_mfd.arg7)
#endif
#if ARG>=8 
	 ,arg8(_mfd.arg8)
#endif
#if ARG>=9 
	 ,arg9(_mfd.arg9)
#endif
#if ARG>=10 
	 ,arg10(_mfd.arg10)
#endif
    {};
    MemFunData& operator=(const MemFunData& _mfd)
    {
      mem_fun_ptr = _mfd.mem_fun_ptr;
      object = _mfd.object;
#if ARG>=1 
      arg1 = _mfd.arg1;
#endif
#if ARG>=2 
      arg2 = _mfd.arg2;
#endif
#if ARG>=3 
      arg3 = _mfd.arg3;
#endif
#if ARG>=4 
      arg4 = _mfd.arg4;
#endif
#if ARG>=5 
      arg5 = _mfd.arg5;
#endif
#if ARG>=6 
      arg6 = _mfd.arg6;
#endif
#if ARG>=7 
      arg7 = _mfd.arg7;
#endif
#if ARG>=8 
      arg8 = _mfd.arg8;
#endif
#if ARG>=9 
      arg9 = _mfd.arg9;
#endif
#if ARG>=10 
      arg10 = _mfd.arg10;
#endif
      return *this;
    };
    static void * thread_entry(void *arg_ptr) {
      MemFunData * mfd = reinterpret_cast<MemFunData *>(arg_ptr);
      ((mfd->object)->*(mfd->mem_fun_ptr))(
#if ARG>=1 
					   mfd->arg1
#endif
#if ARG>=2 
					   , mfd->arg2
#endif
#if ARG>=3 
					   , mfd->arg3
#endif
#if ARG>=4 
					   , mfd->arg4
#endif
#if ARG>=5 
					   , mfd->arg5
#endif
#if ARG>=6 
					   , mfd->arg6
#endif
#if ARG>=7 
					   , mfd->arg7
#endif
#if ARG>=8 
					   , mfd->arg8
#endif
#if ARG>=9 
					   , mfd->arg9
#endif
#if ARG>=10 
					   , mfd->arg10
#endif
					   );
      return NULL;
    };

  };

    class MemArgCollector
  {
  public:
    typedef MemFunData FunDataType;
    MemFunPtr mem_fun_ptr;

    MemArgCollector (MemFunPtr _mem_fun_ptr_) {
      mem_fun_ptr = _mem_fun_ptr_;
    };
        
    MemFunData collectArgs (Class * object
#if ARG>=1 
	    ,Arg1 _arg1
#endif
#if ARG>=2 
	    ,Arg2 _arg2
#endif
#if ARG>=3 
	    ,Arg3 _arg3
#endif
#if ARG>=4 
	    ,Arg4 _arg4
#endif
#if ARG>=5 
	    ,Arg5 _arg5
#endif
#if ARG>=6 
	    ,Arg6 _arg6
#endif
#if ARG>=7 
	    ,Arg7 _arg7
#endif
#if ARG>=8 
	    ,Arg8 _arg8
#endif
#if ARG>=9 
	    ,Arg9 _arg9
#endif
#if ARG>=10 
	    ,Arg10 _arg10
#endif
	    ) {
      return MemFunData(mem_fun_ptr
			,object
#if ARG>=1 
			, _arg1
#endif
#if ARG>=2 
			, _arg2
#endif
#if ARG>=3 
			, _arg3
#endif
#if ARG>=4 
			, _arg4
#endif
#if ARG>=5 
			, _arg5
#endif
#if ARG>=6 
			, _arg6
#endif
#if ARG>=7 
			, _arg7
#endif
#if ARG>=8 
			, _arg8
#endif
#if ARG>=9 
			, _arg9
#endif
#if ARG>=10 
			, _arg10
#endif
			);
    };
  };
};

template <typename Return
#if ARG>=1
,typename Arg1
#endif
#if ARG>=2
,typename Arg2
#endif
#if ARG>=3
,typename Arg3
#endif
#if ARG>=4
,typename Arg4
#endif
#if ARG>=5
,typename Arg5
#endif
#if ARG>=6
,typename Arg6
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
>
typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
>::ArgCollector encapsulate(Return (*_fun_ptr)(
#if ARG>=1
Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
)) {
    return typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
>::ArgCollector(_fun_ptr);
};

template <typename Return
#if ARG>=1
,typename Arg1
#endif
#if ARG>=2
,typename Arg2
#endif
#if ARG>=3
,typename Arg3
#endif
#if ARG>=4
,typename Arg4
#endif
#if ARG>=5
,typename Arg5
#endif
#if ARG>=6
,typename Arg6
#endif
#if ARG>=7 
,typename Arg7 
#endif
#if ARG>=8 
,typename Arg8
#endif
#if ARG>=9 
,typename Arg9
#endif
#if ARG>=10 
,typename Arg10
#endif
,typename Class>
typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
,Class>::MemArgCollector encapsulate(Return (Class::*_mem_fun_ptr)(
#if ARG>=1
Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
)) {
  return typename CLASSNAME<Return
#if ARG>=1
, Arg1
#endif
#if ARG>=2
, Arg2
#endif
#if ARG>=3
, Arg3
#endif
#if ARG>=4
, Arg4
#endif
#if ARG>=5
, Arg5
#endif
#if ARG>=6
, Arg6
#endif
#if ARG>=7
, Arg7
#endif
#if ARG>=8
, Arg8
#endif
#if ARG>=9
, Arg9
#endif
#if ARG>=10
, Arg10
#endif
,Class>::MemArgCollector(_mem_fun_ptr);
};


/**
 * end of file
 * 
 */
