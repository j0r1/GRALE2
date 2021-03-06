#ifndef OPENCLKERNEL_H

#define OPENCLKERNEL_H

#include <errut/errorbase.h>
#include <string>
#include <stdint.h>

typedef void* cl_command_queue;
typedef uint64_t cl_command_queue_properties;
typedef void* cl_context;
typedef void* cl_context_properties;
typedef void* cl_device_id;
typedef uint64_t cl_device_type;
typedef int32_t cl_int;
typedef void* cl_kernel;
typedef void* cl_platform_id;
typedef void* cl_program;
typedef uint32_t cl_uint;
typedef uint32_t cl_program_build_info;
typedef void* cl_mem;
typedef uint64_t cl_mem_flags;
typedef void* cl_event;
typedef uint32_t cl_bool;
typedef float cl_float;
typedef struct { float x, y; } cl_float2;

#define CL_SUCCESS 0
#define CL_DEVICE_TYPE_GPU (1 << 2)
#define CL_PROGRAM_BUILD_LOG 0x1183
#define CL_MEM_WRITE_ONLY                           (1 << 1)
#define CL_MEM_READ_ONLY                            (1 << 2)
#define CL_MEM_COPY_HOST_PTR                        (1 << 5)

class OpenCLKernel : public errut::ErrorBase
{
public:
	OpenCLKernel();
	~OpenCLKernel();

	bool init(const std::string &libraryName);
	bool loadKernel(const std::string &program, const std::string &kernelName, std::string &failLog);
	bool destroy();

	cl_context getContext()										{ return m_context; }
	cl_kernel getKernel() 										{ return m_kernel; }
	cl_command_queue getCommandQueue()		 						{ return m_queue; }

	cl_int (*clBuildProgram)(cl_program program, cl_uint num_devices, const cl_device_id *device_list, const char *options, void (*pfn_notify)(cl_program, void *user_data), void *user_data);
	cl_command_queue (*clCreateCommandQueue)(cl_context context, cl_device_id device, cl_command_queue_properties properties, cl_int *errcode_ret);
	cl_context (*clCreateContext)(cl_context_properties *properties, cl_uint num_devices, const cl_device_id *devices, void *pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data),	void *user_data, cl_int *errcode_ret);
	cl_kernel (*clCreateKernel)(cl_program  program, const char *kernel_name, cl_int *errcode_ret);
	cl_program (*clCreateProgramWithSource)(cl_context context,	cl_uint count, const char **strings, const size_t *lengths, cl_int *errcode_ret);
	cl_int (*clGetDeviceIDs)(cl_platform_id platform, cl_device_type device_type, cl_uint num_entries, cl_device_id *devices, cl_uint *num_devices);
	cl_int (*clGetPlatformIDs)(cl_uint num_entries, cl_platform_id *platforms, cl_uint *num_platforms);
	cl_int (*clGetProgramBuildInfo)(cl_program program, cl_device_id device, cl_program_build_info param_name, size_t param_value_size, void *param_value, size_t *param_value_size_ret);
	cl_int (*clReleaseCommandQueue)(cl_command_queue command_queue);
	cl_int (*clReleaseContext)(cl_context context);
	cl_int (*clReleaseKernel)(cl_kernel kernel);
	cl_int (*clReleaseProgram)(cl_program program);
	
	cl_mem (*clCreateBuffer)(cl_context context, cl_mem_flags flags, size_t size, void *host_ptr, cl_int *errcode_ret);
	cl_int (*clSetKernelArg)(cl_kernel kernel, cl_uint arg_index, size_t arg_size, const void *arg_value);
	cl_int (*clReleaseMemObject)(cl_mem memobj);
	cl_int (*clEnqueueNDRangeKernel)(cl_command_queue command_queue, cl_kernel kernel, cl_uint work_dim, const size_t *global_work_offset, const size_t *global_work_size, const size_t *local_work_size, cl_uint num_events_in_wait_list, const cl_event *event_wait_list, cl_event *event);
	cl_int (*clFinish)(cl_command_queue command_queue);
	cl_int (*clEnqueueReadBuffer)(cl_command_queue command_queue, cl_mem buffer, cl_bool blocking_read, size_t offset, size_t cb, void *ptr, cl_uint num_events_in_wait_list, const cl_event *event_wait_list, cl_event *event);
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

	void *m_pModule;
};

#endif // OPENCLKERNEL_H
