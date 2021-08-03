#pragma once

#include <errut/errorbase.h>
#include <string>
#include <memory>
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
typedef cl_uint cl_event_info;
typedef uint32_t cl_bool;
typedef float cl_float;
typedef struct { float x, y; } cl_float2;
typedef cl_uint cl_platform_info;
typedef cl_uint cl_device_info;

#define CL_SUCCESS 0
#define CL_DEVICE_TYPE_GPU (1 << 2)
#define CL_PROGRAM_BUILD_LOG 0x1183
#define CL_MEM_READ_WRITE                           (1 << 0)
#define CL_MEM_WRITE_ONLY                           (1 << 1)
#define CL_MEM_READ_ONLY                            (1 << 2)
#define CL_MEM_COPY_HOST_PTR                        (1 << 5)
#define CL_PLATFORM_NAME                            0x0902
#define CL_DEVICE_NAME								0x102B
#define CL_COMPLETE									0
#define CL_EVENT_COMMAND_TYPE                       0x11D1
#define CL_EVENT_COMMAND_EXECUTION_STATUS           0x11D3

#ifdef _WIN32
	#define CL_CALLBACK __stdcall
#else
	#define CL_CALLBACK
#endif

class OpenCLLibrary : public errut::ErrorBase
{
public:
	OpenCLLibrary();
	~OpenCLLibrary();

	static std::string getLibraryName();

	bool loadLibrary(const std::string &libraryName);
	bool isOpen() const { return (m_pModule)?true:false; }
	int getDeviceCount() const;

	cl_int (*clBuildProgram)(cl_program program, cl_uint num_devices, const cl_device_id *device_list, const char *options, void (*pfn_notify)(cl_program, void *user_data), void *user_data);
	cl_command_queue (*clCreateCommandQueue)(cl_context context, cl_device_id device, cl_command_queue_properties properties, cl_int *errcode_ret);
	cl_context (*clCreateContext)(cl_context_properties *properties, cl_uint num_devices, const cl_device_id *devices, void *pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data),	void *user_data, cl_int *errcode_ret);
	cl_kernel (*clCreateKernel)(cl_program  program, const char *kernel_name, cl_int *errcode_ret);
	cl_program (*clCreateProgramWithSource)(cl_context context,	cl_uint count, const char **strings, const size_t *lengths, cl_int *errcode_ret);
	cl_int (*clGetDeviceIDs)(cl_platform_id platform, cl_device_type device_type, cl_uint num_entries, cl_device_id *devices, cl_uint *num_devices);
	cl_int (*clGetPlatformIDs)(cl_uint num_entries, cl_platform_id *platforms, cl_uint *num_platforms);
	cl_int (*clGetPlatformInfo)(cl_platform_id platform, cl_platform_info param_name, size_t param_value_size, void *param_value, size_t *param_value_size_ret);
	cl_int (*clGetDeviceInfo)(cl_device_id device, cl_device_info param_name, size_t param_value_size, void *param_value, size_t *param_value_size_ret);

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
	cl_int (*clFlush)(cl_command_queue command_queue);
	cl_int (*clEnqueueReadBuffer)(cl_command_queue command_queue, cl_mem buffer, cl_bool blocking_read, size_t offset, size_t cb, void *ptr, cl_uint num_events_in_wait_list, const cl_event *event_wait_list, cl_event *event);
	cl_int (*clEnqueueWriteBuffer)(cl_command_queue command_queue, cl_mem buffer, cl_bool blocking_write, size_t offset, size_t cb, const void *ptr, cl_uint num_events_in_wait_list, const cl_event *event_wait_list, cl_event *event);
	cl_int (*clSetEventCallback)(cl_event event, cl_int command_exec_callback_type, void (CL_CALLBACK * pfn_notify)(cl_event event, cl_int event_command_status, void *user_data), void *user_data);
	cl_int (*clGetEventInfo)(cl_event event, cl_event_info param_name, size_t param_value_size, void *param_value, size_t *param_value_size_ret);
protected:
	bool getPlatformAndDeviceCount(cl_platform_id &platformId, int &deviceCount) const;
	static std::string getCLErrorString(int errNum);
private:
	void *m_pModule;
};
