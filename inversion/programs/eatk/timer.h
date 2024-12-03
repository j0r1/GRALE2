#pragma once

#include <chrono>

class Timer {
public:
	Timer() {
		start();
	}

	void start() {
		beg = std::chrono::steady_clock::now();
	}

	void stop() {
		end = std::chrono::steady_clock::now();
	}

	double duration() {
		return std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count();
	}
private:
	std::chrono::time_point<std::chrono::steady_clock> beg;
	std::chrono::time_point<std::chrono::steady_clock> end;
};

