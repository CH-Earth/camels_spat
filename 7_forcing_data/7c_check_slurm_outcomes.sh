#!/bin/bash

# Main job ID
MAIN_JOB_ID=740572

# Get detailed information for all array jobs
sacct -j ${MAIN_JOB_ID} --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,Elapsed,NodeList > array_jobs_info.txt

# Filter completed and failed jobs
grep -E "COMPLETED|FAILED|CANCELLED" array_jobs_info.txt > filtered_jobs_info.txt

# Summarize nodes for completed and failed jobs
echo "Summary of Nodes for COMPLETED jobs:"
grep "COMPLETED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq -c

echo "Summary of Nodes for FAILED jobs:"
grep "FAILED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq -c

echo "Summary of Nodes for CANCELLED jobs:"
grep "CANCELLED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq -c

# Check if all COMPLETED jobs occurred on specific nodes
COMPLETED_NODES=$(grep "COMPLETED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq)
FAILED_NODES=$(grep "FAILED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq)
CANCELLED_NODES=$(grep "CANCELLED" filtered_jobs_info.txt | awk '{print $NF}' | sort | uniq)

echo
echo "Nodes where COMPLETED jobs ran:"
echo "${COMPLETED_NODES}"
echo
echo "Nodes where FAILED jobs ran:"
echo "${FAILED_NODES}"
echo
echo "Nodes where CANCELLED jobs ran:"
echo "${CANCELLED_NODES}"
echo

# Determine if all COMPLETED jobs occur on specific nodes
if [ $(echo "${COMPLETED_NODES}" | wc -l) -eq 1 ]; then
    echo "All COMPLETED jobs occurred on the same node(s):"
    echo "${COMPLETED_NODES}"
else
    echo "COMPLETED jobs occurred on multiple nodes."
fi

# Determine if all FAILED jobs occur on specific nodes
if [ $(echo "${FAILED_NODES}" | wc -l) -eq 1 ]; then
    echo "All FAILED jobs occurred on the same node(s):"
    echo "${FAILED_NODES}"
else
    echo "FAILED jobs occurred on multiple nodes."
fi

# Determine if all CANCELLED jobs occur on specific nodes
if [ $(echo "${CANCELLED_NODES}" | wc -l) -eq 1 ]; then
    echo "All CANCELLED jobs occurred on the same node(s):"
    echo "${CANCELLED_NODES}"
else
    echo "CANCELLED jobs occurred on multiple nodes."
fi
