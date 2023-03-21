#pragma once
#include<ostream>
template<typename T> class basic_stopwatch : T {
	typedef typename T BaseTimer;
public:
	// create, optionally start timming an activity
	explicit basic_stopwatch(bool start):
		basic_stopwatch("Stopwatch", start)
	{
	}

	explicit basic_stopwatch(char const* activity = "Stopwatch",
							 bool start = true) :
		m_activity(activity),
		m_log(std::cout)
	{
		if (start)
		{
			Start();
			m_log << "Timer Start:" << BaseTimer::GetMs() << "ms\n";
		}
	}

	// stop and destroy a stop watch
	~basic_stopwatch() {
		m_log << "Timer Stop:" << BaseTimer::GetMs() << "ms\n\n";
	}

	// get last lap time (time of last loop)
	unsigned LapGet() const {
		return BaseTimer::GetMs();
	}

	// predicate: return true if the stopwatch is running
	bool IsStarted() const {
		return BaseTimer::IsStarted();
	}

	// show accumulated time, keep running, set/return lap
	unsigned Show(char const* event = "show") {
		m_activity = event;
		m_log << "Timer Now:" << BaseTimer::GetMs() << "ms\n";
		return BaseTimer::GetMs();
	}

	// (re)start a stopwatch, set/return lap time
	unsigned Start(char const* event = "start") {
		m_activity = event;
		BaseTimer::Start();	
		return BaseTimer::GetMs();
	}

	// stop a running stopwatch, set/return lap time
	unsigned Stop(char const* event = "stop") {
		BaseTimer::Stop();
		m_lap = BaseTimer::GetMs();
		return BaseTimer::GetMs();
	}

private: //Members
	char const*		m_activity;	// "activity" string
	unsigned		m_lap;		// lap time (time of last stop)
	std::ostream&	m_log;		// stream on which to log events	
};