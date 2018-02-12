#ifndef OPENCLKERNEL_H

#define OPENCLKERNEL_H

#include <errut/errorbase.h>
#ifdef GRALE_OPENCL_OPENCL_CL_H
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif // GRALE_OPENCL_OPENCL_CL_H
#include <string>

class OpenCLKernel : public errut::ErrorBase
{
public:
	OpenCLKernel();
	~OpenCLKernel();

	bool init();
	bool loadKernel(const std::string &program, const std::string &kernelName, std::string &failLog);
	bool destroy();

	cl_context getContext()										{ return m_context; }
	cl_kernel getKernel() 										{ return m_kernel; }
	cl_command_queue getCommandQueue()		 						{ return m_queue; }
private:
	static std::string getCLErrorString(int errNum);
	void releaseAll();

	bool m_init;
	std::string m_currentProgram, m_currentKernel;
	
	cl_context m_context;
	cl_command_queue m_queue;
	cl_program m_program;
	cl_kernel m_kernel;
	cl_device_id *m_pDevices;
	int m_deviceIndex;
};

#endif // OPENCLKERNEL_H
