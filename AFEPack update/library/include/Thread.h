/**
 * @file   Thread.h
 * @author Ruo Li
 * @date   Sun Feb 13 21:57:00 2005
 * 
 * @brief  for convenient using of POSIX thread in C++
 * 
 * 
 */

#ifdef MULTITHREAD

#ifndef __Thread_h__
#define __Thread_h__

#include <pthread.h>
#include <iostream>
#include <list>

struct NullThread {};

#define ARG 0
#define CLASSNAME Thread0
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 1
#define CLASSNAME Thread1
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 2
#define CLASSNAME Thread2
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 3
#define CLASSNAME Thread3
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 4
#define CLASSNAME Thread4
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 5
#define CLASSNAME Thread5
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 6
#define CLASSNAME Thread6
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 7
#define CLASSNAME Thread7
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 8
#define CLASSNAME Thread8
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 9
#define CLASSNAME Thread9
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

#define ARG 10
#define CLASSNAME Thread10
#include "ThreadBase.h"
#undef CLASSNAME
#undef ARG

/** 
 * �������ڵ��̸߳�����������1ʱ�������Щ���ֻ��Զ����ж�
 * �̻߳���
 * 
 * @param int �̸߳���
 */
void setThread(unsigned int);
/** 
 * ȡ�����ڵ��̸߳���
 * 
 * 
 * @return �̸߳���
 */
unsigned int getThread();

/*!
\class ThreadManager
\brief �����̹߳���

�����̵߳Ĵ����Լ�����Ĺ�������ʵ����ֻ��д��һ��join������

*/
class ThreadManager
{
 private:
  std::list<pthread_t> threads;	/**< �̵߳ľ�� */
  std::list<void *>    args;	/**< argments of the thread */
  bool                 is_joined; /**< �Ƿ��Ѿ�join���� */
 public:
  ThreadManager() : 
    is_joined(false) {};
  ~ThreadManager() {
    if (!is_joined && !threads.empty()) {
      std::cerr << "Thread manager is not joined before destory."
		<< std::endl;
      abort();
    }
  };
  /** 
   * ����һ���߳�
   * 
   * @param fun_data ���ݸ��߳���ڳ���Ĳ�����
   * @param attr �߳�����
   * 
   * @return �ɹ�����0������-1
   */
  template <class T>
  int start(T fun_data,
	    pthread_attr_t * attr = NULL) {
    pthread_t th;
    T * fun_data_copy = new T(fun_data); /** 
					  * it looks that it is more safe 
					  * if we make a copy here.
					  */
    int error_number = pthread_create(&th, 
				      attr,
				      fun_data.thread_entry, 
				      (void *)fun_data_copy);
    if (error_number == 0) {
      threads.push_back(th);
      args.push_back(fun_data_copy);
#ifdef DEBUG
      //std::cerr << "new thread created successfully. " 
      //	<< std::endl;
#endif
      return 0;
    }
    else {
      std::cout << "thread creating failure with error_number "
		<< error_number << std::endl;
      exit(-1);
    }
  };
  /** 
   * �ȴ��߳�ִ�н���
   * 
   * 
   * @return �ɹ�����0������-1
   */
  template <class T>
  int join(T arg_collector) {
    typedef typename T::FunDataType FunData;
    std::list<pthread_t>::iterator the_th = threads.begin();
    std::list<pthread_t>::iterator end_th = threads.end();
    std::list<void *>::iterator the_arg = args.begin();
    for (;the_th != end_th;the_th ++, the_arg ++) {
      int error_number = pthread_join(*the_th, NULL);
      if (error_number) {
	std::cout << "thread join error with error_number "
		  << error_number << std::endl;
	exit(-1);
      }
      delete (FunData *)(*the_arg);
    }
    clean();
    return 0;
  };
 private:
  void clean() {
    threads.clear();
    args.clear();
    is_joined = false;
  };
};

#endif // __Thread_h__

#endif 

/**
 * end of file
 * 
 */
