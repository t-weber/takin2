/**
 * simple stopwatch for debugging
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date jan-2015 - 2018
 * @license GPLv2 or GPLv3
 */

#ifndef __STOPWATCH_H__
#define __STOPWATCH_H__


#include <chrono>
#include <ctime>
#include <cmath>
#include "../string/string.h"


namespace tl {

template<class T = double>
class Stopwatch
{
	public:
		typedef std::chrono::system_clock::time_point t_tp_sys;
		typedef std::chrono::steady_clock::time_point t_tp_st;
		typedef std::chrono::duration<T> t_dur;
		typedef std::chrono::system_clock::duration t_dur_sys;

	protected:
		t_tp_sys m_timeStart;
		t_tp_st m_timeStart_st, m_timeStop_st;

		t_dur m_dur;
		t_dur_sys m_dur_sys;

		T m_dDur = T(0);

	public:
		Stopwatch() = default;
		~Stopwatch() = default;

		void start()
		{
			m_timeStart = std::chrono::system_clock::now();
			m_timeStart_st = std::chrono::steady_clock::now();
		}

		void stop()
		{
			m_timeStop_st = std::chrono::steady_clock::now();

			m_dur = std::chrono::duration_cast<t_dur>(m_timeStop_st-m_timeStart_st);
			m_dur_sys = std::chrono::duration_cast<t_dur_sys>(m_dur);

			m_dDur = T(t_dur::period::num)/T(t_dur::period::den) * T(m_dur.count());
		}

		T GetDur() const
		{
			return m_dDur;
		}

		static std::string to_str(const t_tp_sys& t)
		{
			std::time_t tStart = std::chrono::system_clock::to_time_t(t);
			std::tm tmStart = *std::localtime(&tStart);

			char cTime[256];
			std::strftime(cTime, sizeof cTime, "%a %Y-%b-%d %H:%M:%S %Z", &tmStart);
			return std::string(cTime);
		}

		std::string GetStartTimeStr() const { return to_str(m_timeStart); }
		std::string GetStopTimeStr() const { return to_str(m_timeStart+m_dur_sys); }

		t_tp_sys GetEstStopTime(T dProg) const
		{
			t_tp_st timeStop_st = std::chrono::steady_clock::now();
			t_dur dur = std::chrono::duration_cast<t_dur>(timeStop_st - m_timeStart_st);
			dur *= (T(1)/dProg);

			t_dur_sys dur_sys = std::chrono::duration_cast<t_dur_sys>(dur);
			t_tp_sys tpEnd = m_timeStart + dur_sys;
			return tpEnd;
		}

		std::string GetEstStopTimeStr(T dProg) const
		{
			return to_str(GetEstStopTime(dProg));
		}
};


template<class t_real = double>
std::string get_duration_str_secs(t_real dDur)
{
	int iAgeMS = int((dDur - std::trunc(dDur)) * t_real(1000.));
	int iAge[] = { int(dDur), 0, 0, 0 };	// s, m, h, d
	const int iConv[] = { 60, 60, 24 };
	const char* pcUnit[] = { "s ", "m ", "h ", "d " };

	for(std::size_t i=0; i<sizeof(iAge)/sizeof(iAge[0])-1; ++i)
	{
		if(iAge[i] > iConv[i])
		{
			iAge[i+1] = iAge[i] / iConv[i];
			iAge[i] = iAge[i] % iConv[i];
		}
	}

	bool bHadPrev = 0;
	std::string strAge;
	for(std::ptrdiff_t i=sizeof(iAge)/sizeof(iAge[0])-1; i>=0; --i)
	{
		if(iAge[i] || bHadPrev)
		{
			strAge += tl::var_to_str(iAge[i]) + pcUnit[i];
			bHadPrev = 1;
		}
	}
	/*if(iAgeMS)*/ { strAge += tl::var_to_str(iAgeMS) + "ms "; bHadPrev = 1; }

	return strAge;
}


template<class t_real = double>
std::string get_duration_str(const std::chrono::duration<t_real>& dur)
{
	using t_dur = std::chrono::duration<t_real>;

	t_real dDurSecs = t_real(t_dur::period::num)/t_real(t_dur::period::den)
		* t_real(dur.count());
	return get_duration_str_secs(dDurSecs);
}

}
#endif
