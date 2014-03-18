#! /bin/bash


#job_tracker=`cat MY_NODE_LIST | sort -u | head -2 | tail -1`
job_tracker=`cat MY_NODE_LIST | sort -u | head -1`
echo 'job_tracker = ' $job_tracker
ssh shivaswa@$job_tracker 'bash -s' < stop_job_tracker.sh
name_node=`cat MY_NODE_LIST | sort -u | head -1`
echo 'name_node = ' $name_node
ssh shivaswa@$name_node 'bash -s' < stop_name_node.sh
