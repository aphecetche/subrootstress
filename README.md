# subrootstress
Benchmark for multicore machines based on rootstress

Requirements

- cmake >= 2.18.12

- g++ >= 4.4

- root >= 5.34.30 < 6 (https://root.cern.ch/content/release-53432)

- git (to retrieve subrootstress source)

Installation

	git clone https://github.com/aphecetche/subrootstress.git
	cd subrootstress
    source /path_to_root_installation/bin/thisroot.sh # set env. for root
	make

Execution

	cd subrootstress
	./compute.sh

Report

Results are in the SubRootMarksVsNProc.pdf file

The bench is a mix of CPU (including floating point operations) and I/O stress.
Of interest in the output figure is the value of the y-axis for 1 core (the higher the better) and the amount of performance decrease when going to more cores.

