#!/bin/bash
script_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$script_dir" || exit 1

logfile="$script_dir/1105feff_batch.log"
echo "Batch started at $(date)" > "$logfile"

max_jobs=32
count=0

while IFS= read -r inp; do
    d="$(dirname "$inp")"
    rel="${d#$script_dir/}"
    runlog="${script_dir}/$(echo "$rel" | tr '/' '_')_run.log"

    (
        cd "$d" || {
            {
                flock -x 200
                echo "$(date) [ERROR] Cannot enter $d" >> "$logfile"
            } 200>>"$logfile"
            exit 1
        }

        {
            flock -x 200
            echo ">>> $(date) [START] $d" >> "$logfile"
        } 200>>"$logfile"

        bash ~/JFEFF_FEFF10/feff10.sh > "$runlog" 2>&1
        status=$?

        {
            flock -x 200
            if [ $status -eq 0 ]; then
                echo "$(date) [DONE] $d" >> "$logfile"
            else
                echo "$(date) [FAILED] $d (exit code $status)" >> "$logfile"
            fi
        } 200>>"$logfile"
    ) &

    ((count++))
    if ((count % max_jobs == 0)); then
        wait
    fi
done < <(find "$script_dir" -type f -name "feff.inp")

wait
sync
echo "All FEFF jobs finished at $(date)" >> "$logfile"