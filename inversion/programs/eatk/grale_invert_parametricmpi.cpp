#include <mpi.h>
#include "freeforminversioncommunicator.h"
#include "mpicommunicatortemplate.h"

int main(int argc, char *argv[])
{
	return main_template<FreeFormInversionCommunicator>(argc, argv);
}
