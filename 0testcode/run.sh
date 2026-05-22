#!/bin/bash
#SBATCH -J front_-2  #每次提交的任务名---记得修改哦           

#SBATCH --cpus-per-task=8   # 每个进程占用一个 cpu 核心
#SBATCH -t 5-00:00:00
#SBATCH --gres=gpu:1         # 如果是gpu任务需要在此行定义gpu数量,此处为1
#TIMESTAMP=$(date +%Y%m%d_%H%M%S)  # 定义日志文件名，包括Job ID和时间戳
#SBATCH --output=Alogs/%j-%x.log

# print script contents
echo "------ Script Content ------"
cat $0
echo
echo "----------------------------"

cd front_-2
date
#matlab -nodesktop -nosplash -nodisplay -r RUN_UBASIS_SCRIPT
matlab -nodesktop -nosplash -nodisplay -r UISNR_batch
# matlab -nodesktop -nosplash -nodisplay -r save_data_binary
# matlab -nodisplay -nosplash -nodesktop -r computerSNR
date 