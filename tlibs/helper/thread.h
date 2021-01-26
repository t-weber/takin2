/**
 * Thread helpers
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date aug-2015
 * @license GPLv2 or GPLv3
 */

#ifndef __TLIBS_THREAD_H__
#define __TLIBS_THREAD_H__


//#define USE_OWN_THREADPOOL

#ifndef USE_OWN_THREADPOOL
	#include <boost/asio.hpp>
#endif

#include <future>
#include <thread>
#include <mutex>
#include <list>
#include <functional>
#include <algorithm>
#include <type_traits>
#include <memory>


namespace tl {


/**
 * thread pool
 */
template<class t_func, class t_startfunc = void(void)>
class ThreadPool
{
	public:
		using t_ret = typename std::result_of<t_func&()>::type;
		using t_fut = std::list<std::future<t_ret>>;
		using t_task = std::list<std::packaged_task<t_ret()>>;


	protected:
#ifdef USE_OWN_THREADPOOL
		std::list<std::unique_ptr<std::thread>> m_lstThreads;
		unsigned int m_tp;	// dummy variable to avoid some ifdefs

		// signal to start jobs
		std::promise<void> m_signalStartIn;
		std::future<void> m_signalStartOut = std::move(m_signalStartIn.get_future());
#else
		boost::asio::thread_pool m_tp;
#endif

		std::mutex m_mtx, m_mtxStart;

		// list of wrapped function to be executed
		t_task m_lstTasks;

		// futures with function return values
		t_fut m_lstFutures;

		// function to run before each thread (not task)
		t_startfunc *m_pThStartFunc = nullptr;


	public:
		ThreadPool(unsigned int iNumThreads = std::thread::hardware_concurrency(),
			t_startfunc* pThStartFunc = nullptr)
			: m_tp{iNumThreads}, m_pThStartFunc{pThStartFunc}
		{
#ifdef USE_OWN_THREADPOOL
			// start 'iNumThreads' threads
			for(unsigned int iThread=0; iThread<iNumThreads; ++iThread)
			{
				m_lstThreads.emplace_back(
					std::unique_ptr<std::thread>(new std::thread([this, pThStartFunc, iThread]()
					{
						// callback to invoke before starting job thread
						if(pThStartFunc) (*pThStartFunc)();
						m_signalStartOut.wait();

						while(1)
						{
							std::unique_lock<std::mutex> lock0(m_mtx);

							// is a task available
							if(m_lstTasks.size() > 0)
							{
								// pop task from list
								std::packaged_task<t_ret()> task =
								std::move(m_lstTasks.front());
								m_lstTasks.pop_front();

								lock0.unlock();

								// run start function and task
								CallStartFunc();
								task();
							}
							else
								break;
						}
					})));
			}
#endif
		}


		virtual ~ThreadPool()
		{
			Join();
		}


		/**
		 * add a function to be executed, giving a packaged task and a future.
		 */
		void AddTask(const std::function<t_ret()>& fkt)
		{
			std::packaged_task<t_ret()> task(fkt);
			std::future<t_ret> fut = task.get_future();

			std::lock_guard<std::mutex> lock(m_mtx);
			m_lstTasks.emplace_back(std::move(task));
			m_lstFutures.emplace_back(std::move(fut));

#ifndef USE_OWN_THREADPOOL
			std::packaged_task<t_ret()>* thetask = &m_lstTasks.back();

			boost::asio::post(m_tp, [this, thetask]() -> void
			{
				CallStartFunc();
				(*thetask)();
			});
#endif
		}


		/**
		 * start tasks (does nothing when using boost threadpool)
		 */
		void Start()
		{
#ifdef USE_OWN_THREADPOOL
			m_signalStartIn.set_value();
#endif
		}


		/**
		 * wait for all tasks to be finished
		 */
		void Join()
		{
#ifdef USE_OWN_THREADPOOL
			std::for_each(m_lstThreads.begin(), m_lstThreads.end(),
				[](std::unique_ptr<std::thread>& pThread)
			{
				if(pThread)
					pThread->join();
			});
#else
			m_tp.join();
#endif
		}


		t_fut& GetResults() { return m_lstFutures; }

		t_task& GetTasks() { return m_lstTasks; }


	protected:
		void CallStartFunc()
		{
			// ensure that this is only called per-thread, not per-task
			std::lock_guard<std::mutex> lockStart(m_mtxStart);

			thread_local bool bThreadAlreadySeen{0};
			if(m_pThStartFunc && !bThreadAlreadySeen)
			{
				bThreadAlreadySeen = 1;
				(*m_pThStartFunc)();
			}
		}
};


}
#endif
