#!/usr/bin/env python3
import re
import subprocess as sp
import shlex
import sys
import time
import logging

logger = logging.getLogger("__name__")
logger.setLevel(40)

STATUS_ATTEMPTS = 20

jobid = int(sys.argv[1])
job_status = "running"

# WARNING this currently has no support for task array jobs

for i in range(STATUS_ATTEMPTS):
    # first try qstat to see if job is running
    # we can use `qstat -s pr -u "*"` to check for all running and pending jobs
    try:
        qstat_res = sp.check_output(shlex.split(f"qstat -s pr")).decode().strip()

        # skip the header using [2:]
        res = {
            int(x.split()[0]) : x.split()[4] for x in qstat_res.splitlines()[2:]
        }

        # job is in an unspecified error state
        if "E" in res[jobid]:
            job_status = "failed"
            break

        job_status = "running"
        break

    except sp.CalledProcessError as e:
        logger.error("qstat process error")
        logger.error(e)
    except KeyError as e:
        # if the job has finished it won't appear in qstat and we should check qacct
        # this will also provide the exit status (0 on success, 128 + exit_status on fail)
        try:
            qacct_res = sp.check_output(
                shlex.split(f"qacct -j {jobid}"), stderr=sp.STDOUT
            )

            m = re.search(r"exit_status\s+([0-9]+)", qacct_res.decode())
            if m is None:
                # qacct returned but accounting hasn't caught up yet
                # (no exit_status line). Retry rather than crash.
                if i >= STATUS_ATTEMPTS - 1:
                    job_status = "failed"
                    break
                time.sleep(5)
                continue

            exit_code = int(m.group(1))
            job_status = "success" if exit_code == 0 else "failed"
            break

        except sp.CalledProcessError as e:
            logger.warning("qacct process error")
            logger.warning(e)
            if i >= STATUS_ATTEMPTS - 1:
                job_status = "failed"
                break
            else:
                # qacct can be quite slow to update on large servers
                time.sleep(5)
        pass

print(job_status)
