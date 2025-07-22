#!/bin/bash

# Script per attivare/disattivare il pulser NIM
ADDRESS="aimtti-tgp3152-00"

if [ "$1" == "--on" ]; then
    /eu/aimtti/aimtti-cmd.py --address "$ADDRESS" --list /eu/aimtti-tgp3152/configs/NIM_trg.config
elif [ "$1" == "--off" ]; then
    /eu/aimtti/aimtti-cmd.py --address "$ADDRESS" --list /eu/aimtti-tgp3152/configs/pulser_OFF.config
else
    echo "Uso: $0 --on | --off"
    exit 1
fi