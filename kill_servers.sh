#!/bin/bash

echo "Killing rwaveserver..."
ps -ef | grep rwaveserver | grep -v grep | awk '{print $2}' | xargs -r kill -9

echo "Killing aimtti-server.py..."
ps -ef | grep aimtti-server.py | grep -v grep | awk '{print $2}' | xargs -r kill -9

echo "Done."
