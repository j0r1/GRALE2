#pragma once

#include "opencllibrary.h"

class OpenCLKernel : public OpenCLLibrary
{
public:
	OpenCLKernel();
	~OpenCLKernel();

	bool init(int devIdx = 0);
	bool loadKernel(const std::string &program, const std::string &kernelName, std::string &failLog);
	void destroy();

	cl_context getContext()										{ return m_context; }
	cl_kernel getKernel() 										{ return m_kernel; }
	cl_command_queue getCommandQueue()		 						{ return m_queue; }
private:
	void releaseAll();

	bool m_init;
	std::string m_currentProgram, m_currentKernel;
	
	cl_context m_context;
	cl_command_queue m_queue;
	cl_program m_program;
	cl_kernel m_kernel;
	cl_device_id m_device;
};
