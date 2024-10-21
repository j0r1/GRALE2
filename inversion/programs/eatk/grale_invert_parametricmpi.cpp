#include <mpi.h>
#include "parametricinversioncommunicator.h"
#include "mpicommunicatortemplate.h"

int main(int argc, char *argv[])
{
	return main_template<ParametricInversionCommunicator>(argc, argv);
}
