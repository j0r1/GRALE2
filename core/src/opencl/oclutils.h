#pragma once

#include "graleconfig.h"
#include "opencllibrary.h"
#include <errut/booltype.h>

namespace grale
{

namespace oclutils
{

// The destructor does not auto-deallocate: need to free memory before
// destroying context etc
class CLMem
{
public:
	void dealloc(OpenCLLibrary &cl)
	{
		if (!m_pMem)
			return;
		
		cl.clReleaseMemObject(m_pMem);

		m_pMem = nullptr;
		m_size = 0;
	}

	errut::bool_t realloc(OpenCLLibrary &cl, cl_context ctx, size_t s) // Only reallocates if more memory is requested
	{
		if (s == 0)
			return "Trying to allocate 0 bytes of GPU memory";

		if (s <= m_size)
			return true;

		cl_int err = 0;
		cl_mem pBuf = cl.clCreateBuffer(ctx, CL_MEM_READ_WRITE, s, nullptr, &err);
		if (pBuf == nullptr || err != CL_SUCCESS)
			return "Can't create buffer of size " + std::to_string(s) + " on GPU: code " + std::to_string(err);

		dealloc(cl);

		m_pMem = pBuf;
		m_size = s;
		return true;
	}

	template<class T> 
	errut::bool_t realloc(OpenCLLibrary &cl, cl_context ctx, const std::vector<T> &buffer)
	{
		return realloc(cl, ctx, buffer.size()*sizeof(T));
	}

	errut::bool_t enqueueWriteBuffer(OpenCLLibrary &cl, cl_command_queue queue, const void *pData, size_t s, bool sync = false)
	{
		if (s > m_size)
			return "Size exceeds GPU buffer size";
		cl_int err = cl.clEnqueueWriteBuffer(queue, m_pMem, sync, 0, s, pData, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueuing write buffer: code " + std::to_string(err);
		return true;
	}

	template<class T>
	errut::bool_t enqueueWriteBuffer(OpenCLLibrary &cl, cl_command_queue queue, const std::vector<T> &data, bool sync = false)
	{
		return enqueueWriteBuffer(cl, queue, data.data(), data.size()*sizeof(T), sync);
	}

	errut::bool_t enqueueReadBuffer(OpenCLLibrary &cl, cl_command_queue queue, void *pData, size_t s, cl_event *pDepEvt, cl_event *pEvt, bool sync = false)
	{
		if (s > m_size)
			return "Reading beyond CPU buffer size";
		size_t num = (pDepEvt)?1:0;

		cl_int err = cl.clEnqueueReadBuffer(queue, m_pMem, sync, 0, s, pData, num, pDepEvt, pEvt);
		if (err != CL_SUCCESS)
			return "Error enqueuing read buffer: code " + std::to_string(err);
		return true;
	}

	template<class T>
	errut::bool_t enqueueReadBuffer(OpenCLLibrary &cl, cl_command_queue queue, std::vector<T> &data, cl_event *pDepEvt, cl_event *pEvt, bool sync = false)
	{
		return enqueueReadBuffer(cl, queue, data.data(), data.size()*sizeof(T), pDepEvt, pEvt, sync);
	}

	cl_mem m_pMem = nullptr;
	size_t m_size = 0;
};

} // end namespace oclutils

} // end namespace
