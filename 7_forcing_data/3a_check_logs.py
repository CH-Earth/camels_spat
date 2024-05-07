import glob
logs = glob.glob('/gpfs/mdiops/globalhome/wmk934/HPC/camels_spat/7_forcing_data/slurm_logs/slurm-1712330_*.out')
num_reruns = 0
for log in logs:
    with open(log, 'r') as file:
        for line in file:
            if 'error' in line.lower() and not 'raise' in line.lower():
                print(f'In file {log}:')
                print(f'   {line}')
                num_reruns += 1
print(f'Total number of rerurns needed = {num_reruns}')