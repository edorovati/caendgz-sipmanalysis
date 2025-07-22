#!/bin/bash
echo "Starting servers..."
# running the server for the pulser and the digitizer
/home/eic/bin/rnohup "/eu/aimtti/aimtti-server.py --address aimtti-tgp3152-00"
/home/eic/bin/rnohup /eu/caen-dt5742b/soft/bin/rwaveserver