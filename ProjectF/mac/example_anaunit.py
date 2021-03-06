import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Number of events to run
#events = 100
events = -1 # all

# Specify output root file name
my_proc.set_ana_output_file( sys.argv[1].rsplit("/",2)[1] + "_FDimAna_" + (str(events) if (events > 0) else "all") + ".root" )

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
my_proc.add_process(fmwk.FDimAna())

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
if events > 0:
    my_proc.run(0,events);
else:
    my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
