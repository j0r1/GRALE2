Debugging
=========

This is mostly for myself.

Starting 'gaserver' with gdb in a tmux session
----------------------------------------------

To debug a rare crash that could be replicated by setting the random number
generators to specific values, I needed to start the 'gaserver' process in
a gdb debugging session. On a PC with a desktop and xterm, a similar approach
as the 'gdb' inverter could be used (starts gdb in an xterm window), but on
the cluster where I was working xterm was not available (and running xterm
over an ssh tunnel is not all that fast). Fortunately tmux was available, 
and could be used in a similar way. I created a custom inverter, that started
gdb using a wrapper script so that some additional modules could be loaded.

Here is the inverter code and the script. Note that a tmux session must
be available, to best to start the python code in its own tmux session.

.. code-block::

    import time
    import socket
    import subprocess
    import grale.privutil as privutil
    
    class DBGLocalCSProcessInverter(inverters.Inverter):
        def __init__(self, numProcesses = 1, feedbackObject = None, serverHelperDebugLevel = 0):
    
            self.procs = [ ] # Do this eary, in case __del__ is called sooner than expected
    
            numHelpers = numProcesses
            if numHelpers < 1 or numHelpers > 256:
                raise InverterException("The number of helper processes should be at least one, and at most 256")
    
            # Obtain a port number to use, let's hope it will still be valid in a few seconds
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.bind(('0.0.0.0', 0))
            bindIp, serverPort = s.getsockname()
            del s
    
            super(DBGLocalCSProcessInverter, self).__init__([ "grale_invert_clientserver", "127.0.0.1", str(serverPort) ], 
                                                          "Local client-server ({} helpers)".format(numHelpers),
                                                          feedbackObject=feedbackObject)
    
    
            n = inversion._getModuleName("general")
            moduleDir = inversion._getModuleDirectory(n)
            
            p = subprocess.Popen([
                "tmux", "new-window",
                "/path/to/gdb.sh gaserver -ex 'set args {} {} {}' -ex run".format(0, serverPort, moduleDir) 
                ])
            self.procs.append(p)
            time.sleep(10) # wait a short while before starting to connect the helpers
    
            for i in range(numHelpers):
                p = subprocess.Popen(["gahelper", str(serverHelperDebugLevel), "127.0.0.1", str(serverPort), moduleDir ])
                self.procs.append(p)
    
            time.sleep(0.5) # Wait a short while so that all helpers are completely detected
        
        def destroy(self):
            """Stops the ``gaserver`` and ``gahelper`` programs."""
            for proc in self.procs:
                try:
                    privutil.terminateProcess(proc, feedbackObject = self.feedback)
                except Exception as e:
                    print("Ignoring exception when terminating gaserver or gahelper: " + str(e))
    
            self.procs = []
    
        def __del__(self):
            self.destroy()
    
.. code-block::

    #!/bin/bash -l

    # Used this on my own machine, needed to activate a conda environment
    # before running gdb
    # source ~/anaconda3-2019.03/bin/activate/grale2
    
    # On the cluster, several modules needed to be loaded, and gdb turned
    # out not to be present on the actual nodes (was available on login
    # node though), so I copied it to my own directory.
    # The seed environment variable needed to be set again, as this script
    # essentially starts a new login session.
    echo "Loading modules"
    module load GSL/2.4-GCCcore-6.4.0
    module load iimpi/2018a
    module load Python/3.6.4-intel-2018a
    source $VSC_DATA/ThinKing/PyGraleIntelDBG/bin/activate
    
    echo "Starting gdb"
    export GRALE_DEBUG_SEED=1583449945
    $VSC_DATA/gdb "$@"
    
    echo "Sleeping"
    
