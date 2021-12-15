#pragma once
///	<description>
/// A simple Timer base to get the pass Time value
/// </description>
/// <aurthor> Kurt </aurthor>
/// 
#include <chrono>

using namespace std::chrono;
class TimerBase
{
public:
	// clears the timer
	TimerBase() : m_start(system_clock::time_point::min()){ }

	// clears the timer
	void Clear() {
		m_start = system_clock::time_point::min();
	}

	// returns true if the timer is running
	bool IsStarted() const {
		return (m_start.time_since_epoch() != system_clock::duration(0));
	}

	// start the timer
	void Start() {
		m_start = system_clock::now();
	}

	// get the number of milliseconds since the timer was started
	unsigned long GetMs() {
		if (IsStarted())
		{
			system_clock::duration diff;
			diff = system_clock::now() - m_start;
			return (unsigned)(duration_cast<milliseconds>(diff).count());
		}
		return 0;
	}

private:
	system_clock::time_point m_start;
};

