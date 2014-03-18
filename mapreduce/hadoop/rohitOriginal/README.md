#####Emulator Hadoop MapReduce

These files come from Rohit's user directoy on CCR

    /ifs/user/shivaswa/my_hadoop

I believe _my\_slurm\_script_ is run first with the map set at _map6.py_ and reduce set as _reduce3.py_

Then _my\_slurm\_script\_for\_final\_reduce_ is run with the map set again as map6.py and the reduce set at _reduce5.py_

Note that the emulator is not run inside Hadoop, it is run as a batch process on the cluster before Hadoop is run.
