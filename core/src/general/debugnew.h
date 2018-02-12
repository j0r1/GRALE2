#ifndef GRALE_DEBUGNEW_H

#define GRALE_DEBUGNEW_H

#ifdef GRALE_USE_DEBUG_NEW

	void *operator new(size_t s, const char filename[],int line);
	void *operator new[](size_t s, const char filename[],int line);
	#define new new(__FILE__,__LINE__)
	
#endif // GRALE_USE_DEBUG_NEW

#endif // GRALE_DEBUGNEW_H

