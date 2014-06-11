for i in $(squeue -u kmarcus2 -o "%i"); do sc $i; done
