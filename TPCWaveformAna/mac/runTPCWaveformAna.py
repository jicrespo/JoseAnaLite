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
for x in xrange( len(sys.argv) - 1 ):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("TPCWaveformAnaBook.root");

# Attach an analysis unit ...
my_proc.add_process(fmwk.TPCWaveformAna())

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run() # Run all events
#my_proc.run(0, 100) # Run 100 events

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
