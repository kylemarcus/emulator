#! /bin/bash

#cat $SLURM_JOB_NODELIST | python node_list.py > MY_NODE_LIST
cp ~/my_hadoop/MY_NODE_LIST /user/shivaswa/hadoop-0.20.1/conf

cd /user/shivaswa/hadoop-0.20.1/conf
#cat $SLURM_JOB_NODELIST | python node_list.py > MY_NODE_LIST
name_node=`cat MY_NODE_LIST | sort -u | head -1`
sed -i 's/\/\/.*:/\/\/'$name_node':/g' core-site.xml
#job_tracker=`cat MY_NODE_LIST | sort -u | head -2 | tail -1`
job_tracker=`cat MY_NODE_LIST | sort -u | head -1`
sed -i 's/>.*:/>'$job_tracker':/g' mapred-site.xml

num_nodes=$1
cat MY_NODE_LIST | sort -u | head -n 1 > masters
cat MY_NODE_LIST | sort -u | tail $((1-num_nodes)) > slaves

cd ~/my_hadoop
echo " name_node = " $name_node
ssh shivaswa@$name_node 'bash -s' < boot_name_node.sh
echo "job_tracker = " $job_tracker
ssh shivaswa@$job_tracker 'bash -s' < boot_job_tracker.sh
