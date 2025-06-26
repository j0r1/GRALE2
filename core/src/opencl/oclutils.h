#pragma once

#include "graleconfig.h"
#include "opencllibrary.h"
#include <errut/booltype.h>

#ifdef GRALE_DEBUG_OPENCLSENTINEL
#include <iostream>
#include <stdexcept>
#endif // GRALE_DEBUG_OPENCLSENTINEL

namespace grale
{

namespace oclutils
{

// The destructor does not auto-deallocate: need to free memory before
// destroying context etc
class CLMem
{
#ifdef GRALE_DEBUG_OPENCLSENTINEL
	constexpr static uint8_t sentinelBytes[] = { 0xde, 0xad, 0xbe, 0xef };
	constexpr static size_t numSentinel = sizeof(sentinelBytes);
#endif // GRALE_DEBUG_OPENCLSENTINEL
public:
	void dealloc(OpenCLLibrary &cl)
	{
		if (!m_pMem)
			return;
		
#ifdef GRALE_DEBUG_OPENCLSENTINEL
		auto queue = cl.getSavedQueue();
		uint8_t buf[numSentinel];
		if (cl.clEnqueueReadBuffer(queue, m_pMem, true, m_size, numSentinel, buf, 0, nullptr, nullptr) != CL_SUCCESS)
		{
			std::cerr << "ERROR! Error reading bytes for sentinel check" << std::endl;
			throw std::runtime_error("ERROR! Error reading bytes for sentinel check");
		}

		for (size_t i = 0 ; i < numSentinel ; i++)
		{
			if (buf[i] != sentinelBytes[i])
			{
				std::cerr << "ERROR! Error checking sentinel" << std::endl;
				throw std::runtime_error("ERROR! Error checking sentinel");
			}
		}
		std::cerr << "DEBUG: CLMem::dealloc sentinel OK" << std::endl;
#endif // GRALE_DEBUG_OPENCLSENTINEL

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

#ifdef GRALE_DEBUG_OPENCLSENTINEL
		size_t extraSize = numSentinel;
#else
		size_t extraSize = 0;
#endif // GRALE_DEBUG_OPENCLSENTINEL

		cl_int err = 0;
		cl_mem pBuf = cl.clCreateBuffer(ctx, CL_MEM_READ_WRITE, s + extraSize, nullptr, &err);
		if (pBuf == nullptr || err != CL_SUCCESS)
			return "Can't create buffer of size " + std::to_string(s) + " on GPU: code " + std::to_string(err);

		dealloc(cl);

		m_pMem = pBuf;
		m_size = s;

#ifdef GRALE_DEBUG_OPENCLSENTINEL
		auto queue = cl.getSavedQueue();
		auto ret = cl.clEnqueueWriteBuffer(queue, m_pMem, true, m_size, numSentinel, sentinelBytes, 0, nullptr, nullptr);
		if (ret != CL_SUCCESS)
		{
			std::cerr << "ERROR! Error writing sentinel bytes: " << ret << std::endl;
			throw std::runtime_error("ERROR! Error writing sentinel bytes");
		}
#endif // GRALE_DEBUG_OPENCLSENTINEL

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

	void swap(CLMem &m)
	{
		cl_mem pTmp = m_pMem;
		size_t tmpSize = m_size;

		m_pMem = m.m_pMem;
		m_size = m.m_size;

		m.m_pMem = pTmp;
		m.m_size = tmpSize;
	}

	cl_mem m_pMem = nullptr;
	size_t m_size = 0;
};

} // end namespace oclutils

} // end namespace
